"""Compare transferred hit annotations with curated GenBank features.

Strand, label spacing, and duplicate identical intervals are ignored.
misc_feature intervals shorter than ``tiny_misc_bp`` are ignored because
added-nucleotide and BbsI marks are not curated consistently.
"""

from __future__ import annotations

import re
from pathlib import Path

from Bio import SeqIO

FeatureKey = tuple[str, str, int, int]


def normalize_label(label: str) -> str:
    return re.sub(r"[\s_]+", " ", label.strip().lower())


def feature_keys(
    path: Path,
    *,
    drop_tiny_misc: bool = True,
    tiny_misc_bp: int = 30,
) -> set[FeatureKey]:
    """Return ``(type, normalized label, start, end)`` using 1-based inclusive coordinates."""
    record = SeqIO.read(str(path), "genbank")
    keys: set[FeatureKey] = set()
    for feat in record.features:
        if feat.type == "source":
            continue
        start = int(feat.location.start) + 1
        end = int(feat.location.end)
        length = end - start + 1
        if drop_tiny_misc and feat.type == "misc_feature" and length < tiny_misc_bp:
            continue
        label = ""
        if feat.qualifiers.get("label"):
            label = str(feat.qualifiers["label"][0])
        keys.add((feat.type, normalize_label(label), start, end))
    return keys


def missing_features(predicted: Path, curated: Path, **kwargs) -> set[FeatureKey]:
    """Curated intervals with no predicted interval of the same type, label, and coordinates."""
    return feature_keys(curated, **kwargs) - feature_keys(predicted, **kwargs)
