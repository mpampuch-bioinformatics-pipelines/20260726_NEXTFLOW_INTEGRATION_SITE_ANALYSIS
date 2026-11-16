#!/usr/bin/env bash
# --flanking so igv-reports 1.17.0 embeds the full extract, not BED±500.
set -euo pipefail

script="$(cd "$(dirname "$0")/.." && pwd)/bin/igv_report_flanking.sh"

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

# Same numbers as the cmc103 combined ptg000010l report: BED 49001-63540
# on an 87175 bp extract. flanking/2 must reach base 1 from 0-based start 49000.
assert_eq "$("$script" 49001 63540 87175)" "98000" "full extract from insertion window"

# Window already spans the sequence: no extra flanking.
assert_eq "$("$script" 1 100 100)" "0" "window is the whole sequence"

# Need more on the right than the left.
assert_eq "$("$script" 10 20 100)" "160" "right flank dominates"

echo "igv_report_flanking_test: ok"
