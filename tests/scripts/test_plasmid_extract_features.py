#!/usr/bin/env python3
"""Lightweight checks for plasmid feature extraction (no Nextflow required)."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "bin" / "plasmid_extract_features.py"
PLASMID = ROOT / "TESTS" / "TEST_DATA" / "PLASMIDS" / "cmc_087.gb"
FEATURE_TYPES = {"CDS", "promoter", "terminator", "misc_feature", "protein_bind"}


def main() -> int:
    if not PLASMID.exists():
        print(f"Missing test plasmid: {PLASMID}", file=sys.stderr)
        return 1

    with tempfile.TemporaryDirectory() as tmp:
        prefix = Path(tmp) / "test_plasmid"
        cmd = [
            sys.executable,
            str(SCRIPT),
            str(PLASMID),
            "--prefix",
            str(prefix),
            "--min-length",
            "8",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            print(proc.stdout)
            print(proc.stderr, file=sys.stderr)
            return proc.returncode

        gb = Path(f"{prefix}.gb")
        gff = Path(f"{prefix}.gff3")
        fa = Path(f"{prefix}.features.fasta")
        assert gb.exists() and gb.stat().st_size > 0
        assert gff.exists() and gff.stat().st_size > 0
        assert fa.exists() and fa.stat().st_size > 0

        gff_types = set()
        lengths_ok = True
        with gff.open() as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                assert len(cols) >= 9
                ftype = cols[2]
                start, end = int(cols[3]), int(cols[4])
                assert ftype in FEATURE_TYPES, f"unexpected type {ftype}"
                assert end - start + 1 >= 8, f"short feature in GFF: {line}"
                gff_types.add(ftype)

        n_fasta = sum(1 for line in fa.open() if line.startswith(">"))
        n_gff = sum(1 for line in gff.open() if line.strip() and not line.startswith("#"))
        assert n_fasta == n_gff, f"FASTA/GFF count mismatch: {n_fasta} vs {n_gff}"
        assert n_fasta > 0
        assert "CDS" in gff_types
        print(f"OK: {n_fasta} features; types={sorted(gff_types)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
