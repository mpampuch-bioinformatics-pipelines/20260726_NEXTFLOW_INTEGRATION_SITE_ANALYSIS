#!/usr/bin/env bash
# Cluster BLAST outfmt 6 hits into loci on the same subject sequence.
#
# Hits are the same locus when they lie on one contig and the gap between them
# is at most --merge-distance bp. That distance should be about the length of
# the inserted construct: features of one insert fall inside it, and a second
# insert farther away stays separate. Different contigs are never merged.
set -euo pipefail

merge_distance=5000
flank=50000

while [[ $# -gt 0 ]]; do
    case "$1" in
        --merge-distance)
            merge_distance="$2"
            shift 2
            ;;
        --flank)
            flank="$2"
            shift 2
            ;;
        --)
            shift
            break
            ;;
        -*)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
        *)
            break
            ;;
    esac
done

blast="${1:?BLAST outfmt 6 file required}"

normalized="$(mktemp)"
trap 'rm -f "$normalized"' EXIT

# subject, start, end, query, original sstart, original send
awk -F '\t' 'NF >= 10 && $1 !~ /^#/ {
    ss = $9 + 0
    se = $10 + 0
    lo = ss < se ? ss : se
    hi = ss > se ? ss : se
    printf "%s\t%d\t%d\t%s\t%d\t%d\n", $2, lo, hi, $1, ss, se
}' "$blast" > "$normalized"

echo -e "locus_id\tsubject\taln_start\taln_end\textract_start\textract_end\tn_members\tmembers"

LC_ALL=C sort -k1,1 -k2,2n -k3,3n "$normalized" | awk -F '\t' -v dist="$merge_distance" -v flank="$flank" '
function flush() {
    if (n == 0) return
    locus++
    ext_start = aln_start - flank
    if (ext_start < 1) ext_start = 1
    ext_end = aln_end + flank
    printf "%d\t%s\t%d\t%d\t%d\t%d\t%d\t%s\n", locus, subject, aln_start, aln_end, ext_start, ext_end, n, members
    n = 0
    members = ""
}
{
    s = $1
    lo = $2 + 0
    hi = $3 + 0
    member = $4 ":" $5 "-" $6
    if (n > 0 && (s != subject || lo > aln_end + dist)) flush()
    if (n == 0) {
        subject = s
        aln_start = lo
        aln_end = hi
        members = member
        n = 1
    } else {
        if (hi > aln_end) aln_end = hi
        members = members ";" member
        n++
    }
}
END { flush() }
'
