#!/usr/bin/env python3
"""Checks for BLAST→hit annotation transfer."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "bin" / "annotate_hit_from_blast.py"


def main() -> int:
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        hit = tmp_path / "hit.fa"
        feats = tmp_path / "features.fa"
        blast = tmp_path / "blast.txt"
        prefix = tmp_path / "annotated"

        # Hit window of 100 bp
        hit.write_text(">ptg000001l:1-100\n" + ("ACGT" * 25) + "\n")
        # Feature matches bases 21-40 of the hit (1-based)
        feats.write_text(
            ">feat_0001 type=CDS label=test_cds\n" + ("ACGT" * 5) + "\n"
        )
        # Perfect hit on forward strand: sstart=21 send=40
        blast.write_text(
            "feat_0001\tptg000001l:1-100\t100.000\t20\t0\t0\t1\t20\t21\t40\t0.0\t40.0\n"
        )

        proc = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--hit-fasta",
                str(hit),
                "--features-fasta",
                str(feats),
                "--blast",
                str(blast),
                "--prefix",
                str(prefix),
                "--min-identity",
                "96",
            ],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            print(proc.stdout)
            print(proc.stderr, file=sys.stderr)
            return proc.returncode

        gff = Path(f"{prefix}.gff3")
        gb = Path(f"{prefix}.gb")
        assert gff.exists() and gb.exists()
        gff_text = gff.read_text()
        assert "CDS" in gff_text
        assert "feat_0001" in gff_text
        assert "21\t40" in gff_text
        assert "LOCUS" in gb.read_text()
        print("OK: annotated CDS 21-40 onto hit")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
