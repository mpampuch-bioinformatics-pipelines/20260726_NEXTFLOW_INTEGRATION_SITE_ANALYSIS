#!/usr/bin/env bash
# Liftoff of genes that overlap the syntenic flanks, from GFF3 and GTF.
set -euo pipefail

script="$(cd "$(dirname "$0")/.." && pwd)/bin/lift_reference_annotations.sh"
module load minimap2/2.31 liftoff/1.6.1

fail() {
    echo "FAIL: $*" >&2
    exit 1
}

assert_eq() {
    local got="$1"
    local want="$2"
    local label="$3"
    [[ "$got" == "$want" ]] || fail "$label: got '$got' want '$want'"
}

dna() {
    local length="$1"
    local seed="$2"
    awk -v n="$length" -v seed="$seed" 'BEGIN {
        srand(seed)
        for (i = 1; i <= n; i++) {
            r = int(rand() * 4)
            printf "%s", substr("ACGT", r + 1, 1)
        }
    }'
}

write_reference() {
    local dir="$1"
    local ref left right insert
    ref="$(dna 4000 11)"
    printf '>ref\n%s\n' "$ref" > "$dir/ref.fa"
    left="$(printf '%s' "$ref" | awk '{ print substr($0, 100, 700) }')"
    right="$(printf '%s' "$ref" | awk '{ print substr($0, 800, 700) }')"
    insert="$(dna 400 17)"
    printf '>hit\n%s%s%s\n' "$left" "$insert" "$right" > "$dir/hit.fa"
}

write_gff() {
    cat > "$1" << 'EOF'
##gff-version 3
ref	synth	gene	200	400	.	+	.	ID=gene_flank;Name=gene_flank
ref	synth	mRNA	200	400	.	+	.	ID=rna_flank;Parent=gene_flank
ref	synth	exon	200	400	.	+	.	ID=exon_flank;Parent=rna_flank
ref	synth	CDS	200	400	.	+	0	ID=cds_flank;Parent=rna_flank
ref	synth	gene	700	1000	.	+	.	ID=gene_span;Name=gene_span
ref	synth	mRNA	700	1000	.	+	.	ID=rna_span;Parent=gene_span
ref	synth	exon	700	1000	.	+	.	ID=exon_span;Parent=rna_span
ref	synth	CDS	700	1000	.	+	0	ID=cds_span;Parent=rna_span
ref	synth	gene	3000	3200	.	+	.	ID=gene_out;Name=gene_out
ref	synth	mRNA	3000	3200	.	+	.	ID=rna_out;Parent=gene_out
ref	synth	exon	3000	3200	.	+	.	ID=exon_out;Parent=rna_out
ref	synth	CDS	3000	3200	.	+	0	ID=cds_out;Parent=rna_out
EOF
}

write_gtf() {
    cat > "$1" << 'EOF'
ref	synth	gene	200	400	.	+	.	gene_id "gene_flank";
ref	synth	transcript	200	400	.	+	.	gene_id "gene_flank"; transcript_id "rna_flank";
ref	synth	exon	200	400	.	+	.	gene_id "gene_flank"; transcript_id "rna_flank";
ref	synth	CDS	200	400	.	+	0	gene_id "gene_flank"; transcript_id "rna_flank";
ref	synth	gene	700	1000	.	+	.	gene_id "gene_span";
ref	synth	transcript	700	1000	.	+	.	gene_id "gene_span"; transcript_id "rna_span";
ref	synth	exon	700	1000	.	+	.	gene_id "gene_span"; transcript_id "rna_span";
ref	synth	CDS	700	1000	.	+	0	gene_id "gene_span"; transcript_id "rna_span";
ref	synth	gene	3000	3200	.	+	.	gene_id "gene_out";
ref	synth	transcript	3000	3200	.	+	.	gene_id "gene_out"; transcript_id "rna_out";
ref	synth	exon	3000	3200	.	+	.	gene_id "gene_out"; transcript_id "rna_out";
ref	synth	CDS	3000	3200	.	+	0	gene_id "gene_out"; transcript_id "rna_out";
EOF
}

status_of() {
    awk -F '\t' -v gene="$2" 'NR > 1 && $1 == gene { print $2 }' "$1"
}

check_synthetic() {
    local label="$1"
    local annotation="$2"
    local lifted="$3"
    local dir="$4"
    "$script" --reference "$dir/ref.fa" --annotation "$annotation" --outdir "$dir/out" --threads 2 "$dir/hit.fa"
    local tsv="$dir/out/hit/validation.tsv"
    [[ -s "$tsv" ]] || fail "$label validation missing"
    [[ -s "$dir/out/hit/$lifted" ]] || fail "$label lifted file missing"
    assert_eq "$(status_of "$tsv" gene_flank)" "lifted" "$label flank gene"
    local span
    span="$(status_of "$tsv" gene_span)"
    [[ "$span" == "partial" || "$span" == "lifted" ]] || fail "$label spanning gene status '$span'"
    if grep -q 'gene_out' "$dir/out/hit/$lifted"; then
        fail "$label lifted the gene outside the extract"
    fi
    awk -F '\t' 'NR > 1 && $1 == "gene_flank" && $4 + 0 >= 0.9 { ok = 1 } END { exit !ok }' "$tsv" \
        || fail "$label flank identity is below 0.9"
}

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

mkdir -p "$tmp/gff" "$tmp/gtf"
write_reference "$tmp/gff"
cp "$tmp/gff/ref.fa" "$tmp/gff/hit.fa" "$tmp/gtf/"
write_gff "$tmp/gff/ann.gff3"
write_gtf "$tmp/gtf/ann.gtf"
gzip -c "$tmp/gff/ann.gff3" > "$tmp/gff/ann.gff3.gz"
gzip -c "$tmp/gtf/ann.gtf" > "$tmp/gtf/ann.gtf.gz"

check_synthetic gff "$tmp/gff/ann.gff3.gz" lifted.gff3 "$tmp/gff"
check_synthetic gtf "$tmp/gtf/ann.gtf.gz" lifted.gtf "$tmp/gtf"

real_hit="/ibex/project/c2303/20260726_NEXTFLOW_INTEGRATION_SITE_ANALYSIS/OUTPUTS/20261007_215728/cmc-087-1-57_C.merolae_unpurged/integration-results/combined/extracted_hits/auto_loci_prediction/1_ptg000011l_441832-542582.fa"
real_ref="/ibex/project/c2303/20260726_NEXTFLOW_INTEGRATION_SITE_ANALYSIS/TESTS/TEST_DATA/REFERENCE-ASSEMBLY/GCF_000091205.1_ASM9120v1_genomic.fna"
ann_dir="/ibex/project/c2303/20260726_NEXTFLOW_INTEGRATION_SITE_ANALYSIS/TEST-DATA/ASSEMBLIES/REFERENCE-ANNOTATION/ASM9120v1"
[[ -s "$real_hit" && -s "$real_ref" ]] || fail "sample 57 hit or reference is missing"

check_real() {
    local label="$1"
    local annotation="$2"
    local lifted_name="$3"
    "$script" --reference "$real_ref" --annotation "$annotation" --outdir "$tmp/real-$label" --threads 4 "$real_hit"
    local hitdir="$tmp/real-$label/1_ptg000011l_441832-542582"
    [[ -s "$hitdir/$lifted_name" ]] || fail "$label lifted file missing"
    [[ -s "$hitdir/validation.tsv" ]] || fail "$label validation missing"
    awk -F '\t' 'NR > 1 && $2 == "lifted" { found = 1 } END { exit !found }' "$hitdir/validation.tsv" \
        || fail "$label lifted no genes"
    local outside outside_id
    outside="$(awk -F '\t' '
        FNR == NR && $1 == "interval" { aligned[$3] = 1; next }
        FNR == NR { next }
        $0 ~ /^#/ || NF < 9 || $3 != "gene" { next }
        !($1 in aligned) { print; exit }
    ' "$hitdir/parsed.tsv" <(gzip -dc "$annotation"))"
    [[ -n "$outside" ]] || fail "$label could not find a gene off the aligned chromosome"
    if [[ "$label" == "gff" ]]; then
        outside_id="$(printf '%s\n' "$outside" | awk -F '\t' '{
            if (match($9, /ID=[^;]+/)) print substr($9, RSTART + 3, RLENGTH - 3)
        }')"
    else
        outside_id="$(printf '%s\n' "$outside" | awk -F '\t' '{
            if (match($9, /gene_id "[^"]+"/)) {
                id = substr($9, RSTART, RLENGTH)
                sub(/^gene_id "/, "", id)
                sub(/"$/, "", id)
                print id
            }
        }')"
    fi
    [[ -n "$outside_id" ]] || fail "$label outside gene has no id"
    if grep -q "$outside_id" "$hitdir/$lifted_name"; then
        fail "$label lifted $outside_id from an unaligned chromosome"
    fi
}

check_real gff "$ann_dir/GCF_000091205.1_ASM9120v1_genomic.gff.gz" lifted.gff3
check_real gtf "$ann_dir/GCF_000091205.1_ASM9120v1_genomic.gtf.gz" lifted.gtf

echo "ok"
