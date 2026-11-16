#!/usr/bin/env bash
# Initial IGV locus on an extracted integration sequence.
#
# Width is max(min_view, insertion + 2*pad), then clamped to 1..length and
# centered on the insertion. Reads can stay filtered to insertion ± pad; this
# window is only the BED that create_report uses as the starting view.
# Prints: view_start view_end
set -euo pipefail

if [[ $# -ne 5 ]]; then
    echo "usage: igv_report_view.sh local_start local_end length min_view pad" >&2
    exit 1
fi

local_start=$1
local_end=$2
len=$3
min_view=$4
pad=$5

for name in local_start local_end len min_view pad; do
    val=${!name}
    if [[ ! "$val" =~ ^[0-9]+$ ]]; then
        echo "${name} must be a non-negative integer, got '${val}'" >&2
        exit 1
    fi
done
if [[ "$local_start" -lt 1 || "$local_end" -lt 1 || "$len" -lt 1 ]]; then
    echo "local_start, local_end, and length must be >= 1" >&2
    exit 1
fi
if [[ "$local_start" -gt "$local_end" ]]; then
    tmp=$local_start
    local_start=$local_end
    local_end=$tmp
fi
if [[ "$local_end" -gt "$len" || "$local_start" -gt "$len" ]]; then
    echo "alignment ${local_start}-${local_end} does not fit sequence length ${len}" >&2
    exit 1
fi

span=$((local_end - local_start + 1))
width=$((span + 2 * pad))
if [[ "$min_view" -gt "$width" ]]; then
    width=$min_view
fi
if [[ "$width" -gt "$len" ]]; then
    width=$len
fi
extra=$((width - span))
if [[ "$extra" -lt 0 ]]; then
    extra=0
fi
left=$((extra / 2))
right=$((extra - left))
view_start=$((local_start - left))
view_end=$((local_end + right))
if [[ "$view_start" -lt 1 ]]; then
    view_end=$((view_end + 1 - view_start))
    view_start=1
fi
if [[ "$view_end" -gt "$len" ]]; then
    view_start=$((view_start - (view_end - len)))
    view_end=$len
fi
if [[ "$view_start" -lt 1 ]]; then
    view_start=1
fi

printf '%d %d\n' "$view_start" "$view_end"
