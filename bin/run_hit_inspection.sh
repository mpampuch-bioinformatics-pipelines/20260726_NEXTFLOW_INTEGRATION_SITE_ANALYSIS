#!/usr/bin/env bash
# Run read inspection, the insertion dot plot, and reference-annotation lift
# into one directory.
#
# The read inspector maps every read it is given and does not subsample.
# This driver subsamples first: the sample 57 HiFi FASTQ is 3.1G, and
# `seqkit sample -n` would load that file into memory. A proportion piped
# to `seqkit head` stops after the requested number of reads.
#
# Workflows, all under bin/ next to this script:
#   read_inspector.sh
#   hit_reference_dotplot.sh
#   lift_reference_annotations.sh
#   paf_hit_reference.awk
#
# Modules: module load liftoff/1.6.1 minimap2/2.31 samtools/1.16.1 htslib/1.16 mummer/4.0 seqkit
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
outdir="$root/OUTPUTS/hit_inspection"
threads="${THREADS:-4}"
sample_n=5000
sample_proportion=0.1
sample_seed=11

reads="/ibex/project/c2303/ALL-BCL-SEQ-READS/BCLCustomers/lauersk/Revio/r84180_20260721_081739/1_D01/version_01/cmc-087-1-57_C.merolae/cmc-087-1-57_C.merolae.hifi_reads.fastq.gz"
reference="$root/TESTS/TEST_DATA/REFERENCE-ASSEMBLY/GCF_000091205.1_ASM9120v1_genomic.fna"
gff="$root/TEST-DATA/ASSEMBLIES/REFERENCE-ANNOTATION/ASM9120v1/GCF_000091205.1_ASM9120v1_genomic.gff.gz"
gtf="$root/TEST-DATA/ASSEMBLIES/REFERENCE-ANNOTATION/ASM9120v1/GCF_000091205.1_ASM9120v1_genomic.gtf.gz"
hits=()

usage() {
    echo "usage: run_hit_inspection.sh [--outdir DIR] [--threads N] [--sample-n N] [hit.fa ...]" >&2
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir)
            outdir="$2"
            shift 2
            ;;
        --threads)
            threads="$2"
            shift 2
            ;;
        --sample-n)
            sample_n="$2"
            shift 2
            ;;
        --reads)
            reads="$2"
            shift 2
            ;;
        --reference)
            reference="$2"
            shift 2
            ;;
        --gff)
            gff="$2"
            shift 2
            ;;
        --gtf)
            gtf="$2"
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
            hits+=("$1")
            shift
            ;;
    esac
done
hits+=("$@")

if [[ ${#hits[@]} -eq 0 ]]; then
    hits=("$root/OUTPUTS/20261007_215728/cmc-087-1-57_C.merolae_unpurged/integration-results/combined/extracted_hits/auto_loci_prediction/1_ptg000011l_441832-542582.fa")
fi

for path in "$reads" "$reference" "$gff" "$gtf" "${hits[@]}"; do
    if [[ ! -s "$path" ]]; then
        echo "Missing or empty input: $path" >&2
        exit 1
    fi
done

# Liftoff's environment ships an older minimap2. Load it first so the
# modules below are the ones at the front of PATH.
if [[ $(type -t module 2>/dev/null) == function || $(type -t module 2>/dev/null) == file ]]; then
    module load liftoff/1.6.1
    module load minimap2/2.31 samtools/1.16.1 htslib/1.16 mummer/4.0 seqkit
fi

mkdir -p "$outdir/reads" "$outdir/read_inspector" "$outdir/dotplot" "$outdir/annotations/gff" "$outdir/annotations/gtf"

# head closes the pipe once it has enough reads, and sample then exits
# with SIGPIPE. That is the stop, not a failed subsample.
set +o pipefail
seqkit sample -p "$sample_proportion" -s "$sample_seed" "$reads" \
    | seqkit head -n "$sample_n" > "$outdir/reads/subset.fastq"
set -o pipefail
if [[ ! -s "$outdir/reads/subset.fastq" ]]; then
    echo "seqkit wrote no reads to $outdir/reads/subset.fastq" >&2
    exit 1
fi

"$root/bin/read_inspector.sh" \
    --outdir "$outdir/read_inspector" \
    --threads "$threads" \
    "$outdir/reads/subset.fastq" \
    "${hits[@]}"

"$root/bin/hit_reference_dotplot.sh" \
    --reference "$reference" \
    --outdir "$outdir/dotplot" \
    --threads "$threads" \
    "${hits[@]}"

"$root/bin/lift_reference_annotations.sh" \
    --reference "$reference" \
    --annotation "$gff" \
    --outdir "$outdir/annotations/gff" \
    --threads "$threads" \
    "${hits[@]}"

"$root/bin/lift_reference_annotations.sh" \
    --reference "$reference" \
    --annotation "$gtf" \
    --outdir "$outdir/annotations/gtf" \
    --threads "$threads" \
    "${hits[@]}"

{
    printf 'outdir\t%s\n' "$outdir"
    printf 'reads\t%s\n' "$reads"
    printf 'subset\t%s\n' "$outdir/reads/subset.fastq"
    printf 'sample_n\t%s\n' "$sample_n"
    printf 'reference\t%s\n' "$reference"
    printf 'gff\t%s\n' "$gff"
    printf 'gtf\t%s\n' "$gtf"
    printf 'workflow\t%s\n' "$root/bin/read_inspector.sh"
    printf 'workflow\t%s\n' "$root/bin/hit_reference_dotplot.sh"
    printf 'workflow\t%s\n' "$root/bin/lift_reference_annotations.sh"
    printf 'workflow\t%s\n' "$root/bin/paf_hit_reference.awk"
    for hit in "${hits[@]}"; do
        printf 'hit\t%s\n' "$hit"
    done
} > "$outdir/inputs.tsv"

echo "$outdir"
