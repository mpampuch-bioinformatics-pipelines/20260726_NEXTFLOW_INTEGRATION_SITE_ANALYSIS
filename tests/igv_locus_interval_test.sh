#!/usr/bin/env bash
# Window on an extracted integration sequence for the IGV report.
set -euo pipefail

script="$(cd "$(dirname "$0")/.." && pwd)/bin/igv_locus_interval.sh"

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

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

# Contig ptg000012l ends at 516922, so blastdbcmd shortens the requested
# extract end 542673. The alignment still sits on the sequence that was returned.
cat > "$tmp/loci.tsv" << 'TSV'
locus_id	subject	aln_start	aln_end	extract_start	extract_end	n_members	members
1	ptg000012l	491923	492673	441923	542673	1	mVenus:491923-492673
2	ptg000012l	100	200	1	100200	1	decoy:100-200
TSV

truncated="$("$script" "$tmp/loci.tsv" ptg000012l_441923-516922 75000 1000)"
assert_eq "$truncated" "49001 51751 50001 50751" "truncated contig"

# Header end equals extract_end. A longer row with the same start must not win.
cat > "$tmp/exact.tsv" << 'TSV'
locus_id	subject	aln_start	aln_end	extract_start	extract_end	n_members	members
1	scaffold_1	150	250	100	300	1	feature:150-250
2	scaffold_1	160	260	100	400	1	other:160-260
TSV

exact="$("$script" "$tmp/exact.tsv" scaffold_1_100-300 201 10)"
assert_eq "$exact" "41 161 51 151" "untruncated window"

# Pad runs off both ends of a short sequence.
cat > "$tmp/clamp.tsv" << 'TSV'
locus_id	subject	aln_start	aln_end	extract_start	extract_end	n_members	members
1	near	5	20	1	100	1	near:5-20
2	far	90	98	1	100	1	far:90-98
TSV

clamp_start="$("$script" "$tmp/clamp.tsv" near_1-100 100 10)"
assert_eq "$clamp_start" "1 30 5 20" "pad clamped at base 1"

clamp_end="$("$script" "$tmp/clamp.tsv" far_1-100 100 10)"
assert_eq "$clamp_end" "80 100 90 98" "pad clamped at sequence length"

if "$script" "$tmp/loci.tsv" missing_1-10 10 0 >/dev/null 2>"$tmp/missing.err"; then
    fail "missing row should exit non-zero"
fi
grep -q 'no locus row' "$tmp/missing.err" || fail "missing row did not report why"

# Two truncated windows share a start and neither end equals the header.
cat > "$tmp/ambiguous.tsv" << 'TSV'
locus_id	subject	aln_start	aln_end	extract_start	extract_end	n_members	members
1	ptg000012l	491923	492673	441923	542673	1	a:491923-492673
2	ptg000012l	491923	492800	441923	560000	1	b:491923-492800
TSV

if "$script" "$tmp/ambiguous.tsv" ptg000012l_441923-516922 75000 1000 >/dev/null 2>"$tmp/ambiguous.err"; then
    fail "ambiguous rows should exit non-zero"
fi
grep -q 'ambiguous locus rows' "$tmp/ambiguous.err" || fail "ambiguous rows did not report why"

echo "igv_locus_interval_test: ok"
