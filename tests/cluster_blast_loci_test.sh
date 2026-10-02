#!/usr/bin/env bash
# Same-contig clustering of the parent combined BLAST hits.
set -euo pipefail

script="$(cd "$(dirname "$0")/.." && pwd)/bin/cluster_blast_loci.sh"

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

# Parent combined hits. Five features span ptg000012l:491801-496743 (max gap 1172 bp).
# mVenus on atg000025l is a different contig and stays its own locus.
cat > "$tmp/parent.blast" << 'EOF'
intron_splice_junction_partialCRE_nucleoplasminNLS_StrepTagII	ptg000012l	100.000	867	0	0	1	867	494352	495218	0.0	1602
Sv40NLS-partialCre_until_intron_splice_junction	ptg000012l	100.000	258	0	0	1	258	493723	493980	7.70e-135	477
CAT	ptg000012l	100.000	648	0	0	1	648	496096	496743	0.0	1197
mVenus	atg000025l	100.000	751	0	0	1	751	15070	14320	0.0	1387
mVenus	ptg000012l	100.000	751	0	0	1	751	491801	492551	0.0	1387
mVenus	ptg000012l	100.000	286	0	0	751	1036	495192	495477	9.03e-150	529
EOF

default_out="$("$script" --merge-distance 5000 --flank 50000 "$tmp/parent.blast")"
default_rows="$(printf '%s\n' "$default_out" | awk 'NR > 1 && NF { c++ } END { print c+0 }')"
assert_eq "$default_rows" "2" "default distance locus count"

printf '%s\n' "$default_out" | awk -F '\t' 'NR > 1 {
    print $1, $2, $3, $4, $5, $6, $7
}' > "$tmp/default.cols"

assert_eq "$(sed -n '1p' "$tmp/default.cols")" "1 atg000025l 14320 15070 1 65070 1" "alternate locus"
assert_eq "$(sed -n '2p' "$tmp/default.cols")" "2 ptg000012l 491801 496743 441801 546743 5" "primary locus"

# A gap of 0 bp merges only overlaps, so the gapped ptg features split
# and the overlapping intron / partial mVenus pair stays together.
zero_out="$("$script" --merge-distance 0 --flank 50000 "$tmp/parent.blast")"
zero_rows="$(printf '%s\n' "$zero_out" | awk 'NR > 1 && NF { c++ } END { print c+0 }')"
assert_eq "$zero_rows" "5" "zero distance locus count"
printf '%s\n' "$zero_out" | grep -q 'intron_splice_junction_partialCRE_nucleoplasminNLS_StrepTagII:494352-495218;mVenus:495192-495477' \
    || fail "overlapping intron and partial mVenus were not merged at distance 0"

# A second interval farther than the merge distance on the same contig is its own locus.
cat "$tmp/parent.blast" > "$tmp/distant.blast"
printf 'other_feature\tptg000012l\t100.000\t100\t0\t0\t1\t100\t600000\t600100\t0.0\t100\n' >> "$tmp/distant.blast"
distant_rows="$("$script" --merge-distance 5000 --flank 50000 "$tmp/distant.blast" | awk 'NR > 1 && NF { c++ } END { print c+0 }')"
assert_eq "$distant_rows" "3" "distant same-contig hit stays separate"

echo "cluster_blast_loci_test: ok"
