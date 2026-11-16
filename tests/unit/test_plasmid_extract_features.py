"""Plasmid extraction keeps 3'UTRs and qualifier text the annotator needs."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).resolve().parents[2]
EXTRACT = ROOT / "bin" / "plasmid_extract_features.py"


def _run(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if result.returncode != 0:
        raise AssertionError(f"{result.stderr}\n{result.stdout}")
    return result


def test_extracts_three_prime_utr_and_cds_qualifiers(tmp_path: Path) -> None:
    plasmid = tmp_path / "plasmid.gb"
    record = SeqRecord(
        Seq("ATGC" * 20),
        id="plasmid",
        name="plasmid",
        description="test plasmid",
    )
    record.annotations["molecule_type"] = "DNA"
    record.features = [
        SeqFeature(
            FeatureLocation(0, 30),
            type="CDS",
            qualifiers={
                "label": ["mVenus optimized"],
                "product": ["yellow protein"],
                "note": ["codon optimized"],
                "codon_start": ["1"],
            },
        ),
        SeqFeature(
            FeatureLocation(30, 80),
            type="3'UTR",
            qualifiers={"label": ["Beta-tubulin 3'UTR (extended)"]},
        ),
        SeqFeature(
            FeatureLocation(0, 6),
            type="misc_feature",
            qualifiers={"label": ["added nucleotides and BbsI site"]},
        ),
    ]
    SeqIO.write(record, str(plasmid), "genbank")

    prefix = tmp_path / "plasmid_out"
    _run(
        [
            sys.executable,
            str(EXTRACT),
            str(plasmid),
            "--prefix",
            str(prefix),
            "--min-length",
            "8",
        ]
    )

    fasta_types = set()
    for rec in SeqIO.parse(str(prefix.with_suffix(".features.fasta")), "fasta"):
        for part in rec.description.split():
            if part.startswith("type="):
                fasta_types.add(part.split("=", 1)[1])
        assert "added_nucleotides" not in rec.description
    assert fasta_types == {"CDS", "3'UTR"}

    meta = json.loads(prefix.with_suffix(".features.json").read_text(encoding="utf-8"))
    by_type = {info["type"]: info for info in meta.values()}
    assert by_type["CDS"]["label"] == "mVenus optimized"
    assert by_type["CDS"]["product"] == "yellow protein"
    assert by_type["CDS"]["note"] == "codon optimized"
    assert by_type["CDS"]["codon_start"] == "1"
    assert by_type["3'UTR"]["label"] == "Beta-tubulin 3'UTR (extended)"
