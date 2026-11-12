#!/usr/bin/env bash
# Indexed contigs and reads for IGV and JBrowse 2.
set -euo pipefail

script="$(cd "$(dirname "$0")/.." && pwd)/bin/read_inspector.sh"
module load minimap2/2.31 samtools/1.16.1 htslib/1.16 seqkit

fail() {
    echo "FAIL: $*" >&2
    exit 1
}

# Deterministic sequence. Avoid a homopolymer so minimap2 can seed it.
dna() {
    local length="$1"
    local seed="$2"
    awk -v n="$length" -v seed="$seed" 'BEGIN {
        srand(seed)
        for (i = 1; i <= n; i++) {
            r = int(rand() * 4)
            printf "%s", substr("ACGT", r + 1, 1)
        }
        printf "\n"
    }'
}

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

seq="$(dna 1500 7)"
printf '>contig\n%s\n' "$seq" > "$tmp/hit.fa"
awk -v seq="$seq" 'BEGIN {
    starts[1] = 1; starts[2] = 400; starts[3] = 800; starts[4] = 1100
    for (i = 1; i <= 4; i++) {
        read = substr(seq, starts[i], 400)
        qual = ""
        for (q = 1; q <= length(read); q++) qual = qual "I"
        print "@read" i
        print read
        print "+"
        print qual
    }
}' > "$tmp/reads.fastq"

"$script" --outdir "$tmp/viewer" --threads 2 "$tmp/reads.fastq" "$tmp/hit.fa"

for required in contigs.fa.gz contigs.fa.gz.fai contigs.fa.gz.gzi reads.bam reads.bam.bai igv_session.xml config.json; do
    [[ -s "$tmp/viewer/$required" ]] || fail "missing $required"
done
samtools quickcheck "$tmp/viewer/reads.bam" || fail "BAM failed quickcheck"
mapped="$(samtools view -c "$tmp/viewer/reads.bam")"
[[ "$mapped" -gt 0 ]] || fail "no mapped reads"
grep -q 'reads.bam' "$tmp/viewer/igv_session.xml" || fail "IGV session does not name the BAM"
grep -q 'BgzipFastaAdapter' "$tmp/viewer/config.json" || fail "JBrowse config is missing the FASTA adapter"
grep -q 'reads.bam' "$tmp/viewer/config.json" || fail "JBrowse config does not name the BAM"
grep -q 'contigs.fa.gz' "$tmp/viewer/config.json" || fail "JBrowse config does not name the FASTA"

# Full HiFi run for sample 57. No downsampling: the FASTQ is the one in the samplesheet.
real_reads="/ibex/project/c2303/ALL-BCL-SEQ-READS/BCLCustomers/lauersk/Revio/r84180_20260721_081739/1_D01/version_01/cmc-087-1-57_C.merolae/cmc-087-1-57_C.merolae.hifi_reads.fastq.gz"
real_hit="/ibex/project/c2303/20260726_NEXTFLOW_INTEGRATION_SITE_ANALYSIS/OUTPUTS/20261007_215728/cmc-087-1-57_C.merolae_unpurged/integration-results/combined/extracted_hits/auto_loci_prediction/1_ptg000011l_441832-542582.fa"
[[ -s "$real_reads" && -s "$real_hit" ]] || fail "sample 57 reads or extracted hit are missing"

# The FASTQ is 3.1G. seqkit sample -n would load it all into memory, so
# stream a fraction and stop after 5000 reads. head closes the pipe and
# sample exits on SIGPIPE; that is how this command stops early.
set +o pipefail
seqkit sample -p 0.1 -s 11 "$real_reads" | seqkit head -n 5000 > "$tmp/subset.fastq"
set -o pipefail
[[ -s "$tmp/subset.fastq" ]] || fail "seqkit subsample is empty"
"$script" --outdir "$tmp/real" --threads 4 "$tmp/subset.fastq" "$real_hit"
samtools quickcheck "$tmp/real/reads.bam" || fail "real BAM failed quickcheck"
real_mapped="$(samtools view -c "$tmp/real/reads.bam")"
[[ "$real_mapped" -gt 0 ]] || fail "real run mapped no reads"
[[ -s "$tmp/real/contigs.fa.gz.fai" && -s "$tmp/real/contigs.fa.gz.gzi" && -s "$tmp/real/reads.bam.bai" ]] || fail "real indexes missing"

echo "ok"
