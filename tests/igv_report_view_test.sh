#!/usr/bin/env bash
set -euo pipefail
script="$(cd "$(dirname "$0")/.." && pwd)/bin/igv_report_view.sh"
fail() { echo "FAIL: $*" >&2; exit 1; }
assert_eq() { [[ "$1" == "$2" ]] || fail "$3: got '$1' want '$2'"; }

# cmc103 combined ptg000010l: 12.5 kb insertion, 23 kb minimum view
assert_eq "$("$script" 50001 62540 87175 23000 1000)" "44771 67770" "min view around cassette"

# Insertion larger than min view: pad each side
assert_eq "$("$script" 1000 50000 87175 23000 1000)" "1 51001" "large insertion plus pad"

# Extract shorter than min view: whole sequence
assert_eq "$("$script" 10 20 100 23000 1000)" "1 100" "clamp to extract"

echo "igv_report_view_test: ok"
