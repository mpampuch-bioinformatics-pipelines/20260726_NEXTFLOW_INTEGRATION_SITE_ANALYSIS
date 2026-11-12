#!/usr/bin/env bash
# The long chromosome shares 4000 bp with the query and the short chromosome
# shares 500 bp, so the best-hit plot must be the long chromosome. The
# whole-reference plot is still written.
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
script="$root/bin/hit_reference_dotplot.sh"

if [[ $(type -t module 2>/dev/null) == function || $(type -t module 2>/dev/null) == file ]]; then
    module load minimap2/2.31 mummer/4.0 samtools/1.16.1
fi

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

python3 - "$tmp" << 'PY'
import random
import sys

out = sys.argv[1]
rng = random.Random(7)

def dna(n):
    return "".join(rng.choice("ACGT") for _ in range(n))

def write(path, records):
    with open(path, "w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start + 60] + "\n")

long_seq = dna(8000)
short_seq = dna(8000)
query = long_seq[500:4500] + short_seq[100:600]
write(f"{out}/reference.fa", [("chrom_long", long_seq), ("chrom_short", short_seq)])
write(f"{out}/hit.fa", [("hit", query)])
PY

"$script" --reference "$tmp/reference.fa" --outdir "$tmp/out" --threads 2 "$tmp/hit.fa"

hitdir="$tmp/out/hit"
chosen="$(tr -d '[:space:]' < "$hitdir/best_chromosome.txt")"
if [[ "$chosen" != "chrom_long" ]]; then
    echo "expected chrom_long, got $chosen" >&2
    exit 1
fi
if [[ ! -s "$hitdir/dotplot.png" ]]; then
    echo "missing whole-reference dotplot.png" >&2
    exit 1
fi
if [[ ! -s "$hitdir/dotplot.best_chromosome.png" ]]; then
    echo "missing dotplot.best_chromosome.png" >&2
    exit 1
fi
if [[ "$(grep -c '^>' "$hitdir/best_chromosome.fa")" -ne 1 ]] || ! grep -q '^>chrom_long' "$hitdir/best_chromosome.fa"; then
    echo "best_chromosome.fa is not only chrom_long" >&2
    exit 1
fi

echo "ok"
