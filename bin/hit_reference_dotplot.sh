#!/usr/bin/env bash
# Map extracted hits back to the reference and draw the insertion.
#
# Later Nextflow passes auto_loci FASTAs unless downstream_hits is all_hits.
# The publication figure is mummerplot of the hit against the whole reference.
# A second mummerplot is the same hit against the reference contig with the
# most aligned bases (sum of PAF column 11). The query coordinates of the
# insertion come from the minimap2 PAF, because mummerplot does not report
# them and Liftoff uses the same interval. modDotPlot is a second identity
# check, installed with uv into the repo-local .venv.
#
# Modules: module load minimap2/2.31 mummer/4.0 samtools/1.16.1
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
paf_awk="$root/bin/paf_hit_reference.awk"
venv="$root/.venv"
reference=""
outdir=""
threads="${THREADS:-4}"
min_insertion=50
max_target_gap=50

usage() {
    echo "usage: hit_reference_dotplot.sh --reference ref.fa --outdir DIR [--threads N] hit.fa [hit.fa...]" >&2
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --reference)
            reference="$2"
            shift 2
            ;;
        --outdir)
            outdir="$2"
            shift 2
            ;;
        --threads)
            threads="$2"
            shift 2
            ;;
        --min-insertion)
            min_insertion="$2"
            shift 2
            ;;
        --max-target-gap)
            max_target_gap="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            break
            ;;
        -*)
            echo "Unknown option: $1" >&2
            usage
            exit 1
            ;;
        *)
            break
            ;;
    esac
done

if [[ -z "$reference" || -z "$outdir" || $# -lt 1 ]]; then
    usage
    exit 1
fi
if [[ ! -s "$reference" ]]; then
    echo "Reference not found or empty: $reference" >&2
    exit 1
fi

# Pip packages stay in this virtualenv. uv creates the interpreter; nothing
# is installed into a module or with system pip.
if [[ ! -x "$venv/bin/moddotplot" ]]; then
    if ! command -v uv >/dev/null 2>&1; then
        echo "uv is required to install moddotplot" >&2
        exit 1
    fi
    if [[ ! -x "$venv/bin/python" ]]; then
        uv venv "$venv"
    fi
    uv pip install --python "$venv/bin/python" moddotplot
fi

mkdir -p "$outdir"

# This mummer 4 build hardcodes its gnuplot launcher to "false", so mummerplot
# writes the .gp script and does not render it. gnuplot 5.4 also rejects the
# "tiny" keyword in that script's png terminal line.
render_mummerplot() {
    local delta="$1"
    local prefix="$2"
    mummerplot --fat --png --layout -p "$prefix" "$delta"
    sed -i 's/set terminal png tiny/set terminal png/' "${prefix}.gp"
    gnuplot "${prefix}.gp"
}

# Reference contig with the most aligned bases. Column 11 is the alignment
# block length; column 10 is the match count, used only when block lengths tie.
# A remaining tie follows name order so the choice does not depend on PAF order.
best_hit_chromosome() {
    local paf="$1"
    awk -F '\t' '
        NF < 12 || $1 ~ /^#/ { next }
        {
            primary = 1
            for (i = 13; i <= NF; i++) if ($i ~ /^tp:A:/ && $i != "tp:A:P") primary = 0
            if (!primary) next
            block[$6] += $11
            matches[$6] += $10
        }
        END {
            best = ""
            for (target in block) {
                if (best == "" \
                    || block[target] > block[best] \
                    || (block[target] == block[best] && matches[target] > matches[best]) \
                    || (block[target] == block[best] && matches[target] == matches[best] && target < best))
                    best = target
            }
            if (best == "") exit 1
            print best
        }
    ' "$paf"
}

for hit in "$@"; do
    if [[ ! -s "$hit" ]]; then
        echo "Hit FASTA not found or empty: $hit" >&2
        exit 1
    fi
    name="$(basename "$hit")"
    name="${name%.gz}"
    name="${name%.fa}"
    name="${name%.fasta}"
    name="${name%.fna}"
    hitdir="$outdir/$name"
    mkdir -p "$hitdir"

    minimap2 -x asm5 -c --secondary=no -t "$threads" "$reference" "$hit" > "$hitdir/align.paf"
    awk -v min_insertion="$min_insertion" -v max_target_gap="$max_target_gap" \
        -f "$paf_awk" "$hitdir/align.paf" > "$hitdir/parsed.tsv"
    {
        printf 'query\ttarget\tstrand\tquery_start\tquery_end\ttarget_start\ttarget_end\tquery_gap\ttarget_gap\n'
        awk -F '\t' '$1 == "gap" { print $2,$3,$4,$5,$6,$7,$8,$9,$10 }' OFS='\t' "$hitdir/parsed.tsv"
    } > "$hitdir/insertion_gaps.tsv"

    if ! chrom="$(best_hit_chromosome "$hitdir/align.paf")"; then
        echo "No primary minimap2 alignment of $hit to $reference" >&2
        exit 1
    fi
    printf '%s\n' "$chrom" > "$hitdir/best_chromosome.txt"
    samtools faidx "$reference" "$chrom" > "$hitdir/best_chromosome.fa"

    (
        cd "$hitdir"
        nucmer --threads "$threads" -p aln "$reference" "$hit"
        # A delta with no cluster header has no alignment to draw.
        if ! grep -q '^>' aln.delta; then
            echo "No nucmer alignment of $hit to $reference" >&2
            exit 1
        fi
        render_mummerplot aln.delta dotplot

        nucmer --threads "$threads" -p aln.best_chromosome best_chromosome.fa "$hit"
        if ! grep -q '^>' aln.best_chromosome.delta; then
            echo "No nucmer alignment of $hit to $chrom" >&2
            exit 1
        fi
        render_mummerplot aln.best_chromosome.delta dotplot.best_chromosome
    )

    # modDotPlot otherwise plots every sequence against every other sequence.
    # On a chromosome-split reference that is the whole genome, not the hit.
    awk -F '\t' '$1 == "interval" { print $3 "\t" $2 }' "$hitdir/parsed.tsv" \
        | sort -u > "$hitdir/moddotplot.pairs"
    if [[ ! -s "$hitdir/moddotplot.pairs" ]]; then
        echo "No reference interval to pass to modDotPlot for $hit" >&2
        exit 1
    fi
    # --pairs needs a FASTA index beside each file it reads. Index copies so
    # the caller's reference is left untouched.
    cp "$reference" "$hitdir/reference.fa"
    cp "$hit" "$hitdir/query.fa"
    samtools faidx "$hitdir/reference.fa"
    samtools faidx "$hitdir/query.fa"
    # -s takes every name after it. A second -s would replace the first.
    mapfile -t sequence_names < <(awk -F '\t' '{ print $1; print $2 }' "$hitdir/moddotplot.pairs" | sort -u)
    "$venv/bin/moddotplot" static \
        -f "$hitdir/reference.fa" "$hitdir/query.fa" \
        -s "${sequence_names[@]}" \
        --pairs "$hitdir/moddotplot.pairs" \
        --compare-only \
        -o "$hitdir/moddotplot"
done
