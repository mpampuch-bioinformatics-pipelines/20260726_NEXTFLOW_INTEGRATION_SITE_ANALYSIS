"""Curated-vs-predicted comparison ignores strand, label spelling, and tiny misc features."""

from __future__ import annotations

import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from annotation_compare import missing_features


def _write_gb(path: Path, features: list[SeqFeature]) -> None:
    record = SeqRecord(Seq("A" * 80), id="hit", name="hit", description="hit window")
    record.annotations["molecule_type"] = "DNA"
    record.features = features
    SeqIO.write(record, str(path), "genbank")


def test_match_ignores_strand_label_spacing_duplicates_and_tiny_misc(tmp_path: Path) -> None:
    curated = tmp_path / "curated.gb"
    predicted = tmp_path / "predicted.gb"
    cds = SeqFeature(
        FeatureLocation(9, 29),
        type="CDS",
        qualifiers={"label": ["mVenus optimized"]},
    )
    tiny = SeqFeature(
        FeatureLocation(0, 8),
        type="misc_feature",
        qualifiers={"label": ["added nucleotides and BbsI site"]},
    )
    _write_gb(curated, [cds, cds, tiny])
    _write_gb(
        predicted,
        [
            SeqFeature(
                FeatureLocation(9, 29, strand=-1),
                type="CDS",
                qualifiers={"label": ["mVenus_optimized"]},
            )
        ],
    )
    assert missing_features(predicted, curated) == set()


def test_reports_a_missing_tandem_copy(tmp_path: Path) -> None:
    curated = tmp_path / "curated.gb"
    predicted = tmp_path / "predicted.gb"
    first = SeqFeature(
        FeatureLocation(9, 29),
        type="terminator",
        qualifiers={"label": ["NOS terminator"]},
    )
    second = SeqFeature(
        FeatureLocation(49, 69),
        type="terminator",
        qualifiers={"label": ["NOS terminator"]},
    )
    _write_gb(curated, [first, second])
    _write_gb(predicted, [first])
    assert missing_features(predicted, curated) == {
        ("terminator", "nos terminator", 50, 69)
    }
