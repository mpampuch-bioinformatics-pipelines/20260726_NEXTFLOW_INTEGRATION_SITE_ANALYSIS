#!/usr/bin/env bash
# Lift reference-genome annotations onto extracted hit contigs.
#
# Later Nextflow runs this only when an annotation file is set, on the same
# hit set as read mapping: auto_loci unless downstream_hits is all_hits.
# minimap2 finds the syntenic reference interval. Liftoff then maps the genes
# that overlap that interval (Shumate and Salzberg, GBE 2021). A gene that
# sits entirely inside one flank alignment must lift with coverage and
# sequence identity at least 0.5. A gene broken by the insertion is partial
# and does not fail the run.
#
# Modules: module load minimap2/2.31 liftoff/1.6.1
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
paf_awk="$root/bin/paf_hit_reference.awk"
reference=""
annotation=""
outdir=""
threads="${THREADS:-4}"
min_coverage=0.5
min_identity=0.5
min_insertion=50
max_target_gap=50

usage() {
    echo "usage: lift_reference_annotations.sh --reference ref.fa --annotation ref.gff3[.gz] --outdir DIR hit.fa [hit.fa...]" >&2
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --reference)
            reference="$2"
            shift 2
            ;;
        --annotation)
            annotation="$2"
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
        --min-coverage)
            min_coverage="$2"
            shift 2
            ;;
        --min-identity)
            min_identity="$2"
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

if [[ -z "$reference" || -z "$annotation" || -z "$outdir" || $# -lt 1 ]]; then
    usage
    exit 1
fi
if [[ ! -s "$reference" || ! -s "$annotation" ]]; then
    echo "Reference and annotation must exist and be non-empty" >&2
    exit 1
fi

case "$annotation" in
    *.gtf|*.gtf.gz) format=gtf ;;
    *.gff|*.gff3|*.gff.gz|*.gff3.gz) format=gff ;;
    *)
        echo "Annotation must be GFF3 or GTF, optionally gzipped: $annotation" >&2
        exit 1
        ;;
esac

mkdir -p "$outdir"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
if [[ "$annotation" == *.gz ]]; then
    gzip -dc "$annotation" > "$work/annotation"
else
    cp "$annotation" "$work/annotation"
fi

failed=0

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
    # Liftoff writes a minimap2 index beside the target FASTA. Copy the hit
    # so that index is not left next to the caller's file.
    cp "$hit" "$hitdir/target.fa"
    hit="$hitdir/target.fa"

    minimap2 -x asm5 -c --secondary=no -t "$threads" "$reference" "$hit" > "$hitdir/align.paf"
    awk -v min_insertion="$min_insertion" -v max_target_gap="$max_target_gap" \
        -f "$paf_awk" "$hitdir/align.paf" > "$hitdir/parsed.tsv"
    if ! awk -F '\t' '$1 == "interval" { found = 1 } END { exit !found }' "$hitdir/parsed.tsv"; then
        echo "No minimap2 alignment of $hit to the reference" >&2
        exit 1
    fi

    subset="$hitdir/subset.$format"
    awk -F '\t' -v format="$format" '
        function overlaps(seq, s, e,    i) {
            if (!(seq in nint)) return 0
            for (i = 1; i <= nint[seq]; i++) {
                if (s <= iend[seq, i] && e >= istart[seq, i]) return 1
            }
            return 0
        }
        function gff_attr(attr, key,    n, parts, i, kv) {
            n = split(attr, parts, ";")
            for (i = 1; i <= n; i++) {
                split(parts[i], kv, "=")
                if (kv[1] == key) return kv[2]
            }
            return ""
        }
        function gtf_attr(attr, key,    re, val) {
            re = key " \"([^\"]*)\""
            if (match(attr, re)) {
                val = substr(attr, RSTART, RLENGTH)
                sub(/^[^"]*"/, "", val)
                sub(/"$/, "", val)
                return val
            }
            return ""
        }
        FNR == NR {
            if ($1 == "interval") {
                seq = $3
                nint[seq]++
                istart[seq, nint[seq]] = $4 + 0
                iend[seq, nint[seq]] = $5 + 0
            }
            next
        }
        {
            nline++
            lines[nline] = $0
            if ($0 ~ /^#/ || NF < 8) next
            seqs[nline] = $1
            starts[nline] = $4 + 0
            ends[nline] = $5 + 0
            types[nline] = $3
            attrs[nline] = $9
        }
        END {
            if (format == "gtf") {
                for (i = 1; i <= nline; i++) {
                    if (types[i] == "gene" && overlaps(seqs[i], starts[i], ends[i])) {
                        gid = gtf_attr(attrs[i], "gene_id")
                        if (gid != "") keep[gid] = 1
                    }
                }
                for (i = 1; i <= nline; i++) {
                    if (lines[i] ~ /^#/) {
                        print lines[i]
                        continue
                    }
                    gid = gtf_attr(attrs[i], "gene_id")
                    if (gid != "" && (gid in keep)) print lines[i]
                }
            } else {
                for (i = 1; i <= nline; i++) {
                    if (types[i] == "gene" && overlaps(seqs[i], starts[i], ends[i])) {
                        gid = gff_attr(attrs[i], "ID")
                        if (gid != "") selected[gid] = 1
                    }
                }
                changed = 1
                while (changed) {
                    changed = 0
                    for (i = 1; i <= nline; i++) {
                        if (attrs[i] == "") continue
                        id = gff_attr(attrs[i], "ID")
                        parent = gff_attr(attrs[i], "Parent")
                        if (id == "" || parent == "" || (id in selected)) continue
                        npar = split(parent, pars, ",")
                        for (p = 1; p <= npar; p++) {
                            if (pars[p] in selected) {
                                selected[id] = 1
                                changed = 1
                                break
                            }
                        }
                    }
                }
                for (i = 1; i <= nline; i++) {
                    if (lines[i] ~ /^#/) {
                        print lines[i]
                        continue
                    }
                    id = gff_attr(attrs[i], "ID")
                    parent = gff_attr(attrs[i], "Parent")
                    if (id != "" && (id in selected)) {
                        print lines[i]
                        continue
                    }
                    npar = split(parent, pars, ",")
                    for (p = 1; p <= npar; p++) {
                        if (pars[p] in selected) {
                            print lines[i]
                            break
                        }
                    }
                }
            }
        }
    ' "$hitdir/parsed.tsv" "$work/annotation" > "$subset"

    if [[ "$format" == "gtf" ]]; then
        # NCBI gene rows set transcript_id to an empty string. gffutils, which
        # Liftoff uses to load the file, indexes that attribute and exits on
        # the empty value. Real transcript ids are left in place.
        sed -i 's/transcript_id ""; *//' "$subset"
    fi

    if [[ "$format" == "gff" ]]; then
        lifted="$hitdir/lifted.gff3"
    else
        lifted="$hitdir/lifted.gtf"
    fi

    if ! awk -F '\t' '$3 == "gene" { found = 1 } END { exit !found }' "$subset"; then
        echo "No annotated genes overlap the syntenic interval of $name" >&2
        printf 'gene_id\tstatus\tcoverage\tsequence_id\tvalid_orfs\tnote\n' > "$hitdir/validation.tsv"
        : > "$hitdir/unmapped_features.txt"
        continue
    fi

    liftoff \
        -g "$subset" \
        -o "$lifted" \
        -u "$hitdir/unmapped_features.txt" \
        -cds \
        -dir "$hitdir/intermediate" \
        -p "$threads" \
        "$hit" \
        "$reference"

    awk -F '\t' -v format="$format" -v min_cov="$min_coverage" -v min_id="$min_identity" \
        -v failfile="$hitdir/validation.failed" '
        function attr_value(s, key,    re, val) {
            re = key "=([^; \t]+)"
            if (match(s, re)) {
                val = substr(s, RSTART, RLENGTH)
                sub(/^[^=]*=/, "", val)
                gsub(/"/, "", val)
                return val
            }
            re = key " \"([^\"]*)\""
            if (match(s, re)) {
                val = substr(s, RSTART, RLENGTH)
                sub(/^[^"]*"/, "", val)
                sub(/"$/, "", val)
                return val
            }
            return ""
        }
        function gene_of(attr) {
            if (format == "gtf") return attr_value(attr, "gene_id")
            id = attr_value(attr, "ID")
            if (id ~ /^gene/) return id
            return ""
        }
        FNR == NR && FILENAME ~ /parsed.tsv$/ {
            if ($1 == "alignment") {
                na++
                aseq[na] = $3
                astart[na] = $7 + 0
                aend[na] = $8 + 0
            } else if ($1 == "gap") {
                nb++
                bseq[nb] = $3
                bpos[nb] = $7 + 0
            }
            next
        }
        FILENAME ~ /subset/ {
            if ($0 ~ /^#/ || NF < 9 || $3 != "gene") next
            if (format == "gtf") gid = attr_value($9, "gene_id")
            else gid = attr_value($9, "ID")
            if (gid == "") next
            ng++
            gene[ng] = gid
            gseq[ng] = $1
            gstart[ng] = $4 + 0
            gend[ng] = $5 + 0
            next
        }
        FILENAME ~ /lifted/ {
            if ($0 ~ /^#/ || NF < 9) next
            if (format == "gtf") gid = attr_value($9, "gene_id")
            else {
                gid = attr_value($9, "ID")
                parent = attr_value($9, "Parent")
            }
            cov = attr_value($9, "coverage")
            sid = attr_value($9, "sequence_ID")
            vorf = attr_value($9, "valid_ORFs")
            if (gid != "" && ($3 == "gene" || (format == "gtf" && $3 == "gene"))) {
                if (cov != "") lcov[gid] = cov
                if (sid != "") lid[gid] = sid
                seen[gid] = 1
            }
            if (format == "gff" && $3 == "gene") {
                if (cov != "") lcov[gid] = cov
                if (sid != "") lid[gid] = sid
                seen[gid] = 1
            }
            if (vorf != "" && gid != "") vorfs[gid] = vorf
            if (format == "gff" && vorf != "" && parent != "") vorfs[parent] = vorf
            next
        }
        END {
            print "gene_id\tstatus\tcoverage\tsequence_id\tvalid_orfs\tnote"
            for (i = 1; i <= ng; i++) {
                gid = gene[i]
                spans = 0
                for (b = 1; b <= nb; b++) {
                    if (bseq[b] == gseq[i] && gstart[i] <= bpos[b] && gend[i] > bpos[b]) spans = 1
                }
                contained = 0
                for (a = 1; a <= na; a++) {
                    if (aseq[a] == gseq[i] && astart[a] <= gstart[i] && gend[i] <= aend[a]) contained = 1
                }
                cov = (gid in lcov) ? lcov[gid] : ""
                sid = (gid in lid) ? lid[gid] : ""
                vorf = (gid in vorfs) ? vorfs[gid] : ""
                note = ""
                if (spans) {
                    status = "partial"
                    note = "spans insertion"
                } else if (!contained) {
                    status = "partial"
                    note = "extends outside a single flank alignment"
                } else if (!(gid in seen) || cov == "" || sid == "" || cov + 0 < min_cov || sid + 0 < min_id) {
                    status = "unmapped"
                    note = "flank gene below coverage or identity"
                    print gid > failfile
                } else {
                    status = "lifted"
                }
                print gid "\t" status "\t" cov "\t" sid "\t" vorf "\t" note
            }
        }
    ' "$hitdir/parsed.tsv" "$subset" "$lifted" > "$hitdir/validation.tsv"

    if [[ -s "$hitdir/validation.failed" ]]; then
        failed=1
    fi
done

exit "$failed"
