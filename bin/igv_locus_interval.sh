#!/usr/bin/env bash
# Padded alignment window on one extracted integration sequence.
#
# loci.tsv is in genome coordinates. The FASTA name is the blastdbcmd range
# subject:extract_start-extract_end with the colon rewritten to an underscore
# so IGV does not parse it as a locus. blastdbcmd shortens that end when the
# requested range runs past the contig, so the name's end can be smaller than
# extract_end. Match on subject and extract_start; an exact extract_end wins
# when several rows share that start.
#
# Prints: keep_start keep_end local_start local_end
# local_* is the alignment on the extracted sequence (1-based, inclusive).
# keep_* is that interval expanded by pad and clamped to 1..length.
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "usage: igv_locus_interval.sh loci.tsv sequence_name length pad" >&2
    exit 1
fi

loci=$1
header=$2
len=$3
pad=$4

if [[ ! -s "$loci" ]]; then
    echo "locus table not found or empty: ${loci}" >&2
    exit 1
fi
if [[ ! "$len" =~ ^[0-9]+$ ]] || [[ "$len" -lt 1 ]]; then
    echo "sequence length must be a positive integer, got '${len}'" >&2
    exit 1
fi
if [[ ! "$pad" =~ ^[0-9]+$ ]]; then
    echo "pad must be a non-negative integer, got '${pad}'" >&2
    exit 1
fi
# Subject ids can contain underscores. The coordinate suffix is the last
# _<start>-<end>, which is the colon blastdbcmd wrote before that suffix.
if [[ ! "$header" =~ ^(.+)_([0-9]+)-([0-9]+)$ ]]; then
    echo "sequence name '${header}' has no trailing _start-end coordinates" >&2
    exit 1
fi

subject="${BASH_REMATCH[1]}"
header_start="${BASH_REMATCH[2]}"
header_end="${BASH_REMATCH[3]}"

if [[ "$header_end" -lt "$header_start" ]]; then
    echo "sequence name '${header}' has end before start" >&2
    exit 1
fi

span=$((header_end - header_start + 1))
if [[ "$span" -ne "$len" ]]; then
    echo "sequence length ${len} does not match ${header} span ${span}" >&2
    exit 1
fi

awk -F '\t' -v subject="$subject" -v header_start="$header_start" -v header_end="$header_end" -v len="$len" -v pad="$pad" '
NR == 1 && $1 == "locus_id" { next }
NF < 6 || $1 ~ /^#/ { next }
$2 == subject && ($5 + 0) == header_start && ($6 + 0) >= header_end {
    n++
    exact[n] = (($6 + 0) == header_end)
    aln_start[n] = $3 + 0
    aln_end[n] = $4 + 0
    extract_start[n] = $5 + 0
}
END {
    chosen = 0
    n_exact = 0
    for (i = 1; i <= n; i++) {
        if (!exact[i]) continue
        n_exact++
        chosen = i
    }
    if (n_exact > 1) {
        print "ambiguous locus rows for " subject " extract_start " header_start > "/dev/stderr"
        exit 1
    }
    if (n_exact == 0) {
        if (n == 0) {
            print "no locus row for " subject " starting at " header_start " and covering " header_end > "/dev/stderr"
            exit 1
        }
        if (n > 1) {
            print "ambiguous locus rows for " subject " extract_start " header_start > "/dev/stderr"
            exit 1
        }
        chosen = 1
    }

    local_start = aln_start[chosen] - extract_start[chosen] + 1
    local_end = aln_end[chosen] - extract_start[chosen] + 1
    if (local_start > local_end) {
        tmp = local_start
        local_start = local_end
        local_end = tmp
    }
    if (local_end < 1 || local_start > len) {
        print "alignment " local_start "-" local_end " does not overlap sequence length " len > "/dev/stderr"
        exit 1
    }

    keep_start = local_start - pad
    if (keep_start < 1) keep_start = 1
    keep_end = local_end + pad
    if (keep_end > len) keep_end = len
    if (keep_start > keep_end) {
        print "padded window is empty for alignment " local_start "-" local_end > "/dev/stderr"
        exit 1
    }

    printf "%d %d %d %d\n", keep_start, keep_end, local_start, local_end
}
' "$loci"
