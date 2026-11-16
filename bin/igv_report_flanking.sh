#!/usr/bin/env bash
# --flanking value for igv-reports 1.17.0 so the embedded slice is the extract.
#
# create_report uses flanking/2 on each side of the 0-based BED, then always
# stores FASTA as that region ± 500. --flanking 0 therefore still leaves the
# chromosome length equal to the extract while sequence exists only on the
# insertion window ± 500. Expand to 1..length; the initial locus stays the BED.
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "usage: igv_report_flanking.sh keep_start keep_end length" >&2
    exit 1
fi

keep_start=$1
keep_end=$2
len=$3

for name in keep_start keep_end len; do
    val=${!name}
    if [[ ! "$val" =~ ^[0-9]+$ ]]; then
        echo "${name} must be a non-negative integer, got '${val}'" >&2
        exit 1
    fi
done
if [[ "$keep_start" -lt 1 || "$keep_end" -lt 1 || "$len" -lt 1 ]]; then
    echo "keep_start, keep_end, and length must be >= 1" >&2
    exit 1
fi

# 0-based BED start is keep_start-1. flanking/2 must reach coordinate 1.
left=$((keep_start - 1))
right=$((len - keep_end))
if [[ "$right" -lt 0 ]]; then
    right=0
fi
each=$left
if [[ "$right" -gt "$each" ]]; then
    each=$right
fi
echo $((each * 2))
