#!/usr/bin/env python3
"""Extract biologically relevant plasmid features to GFF3, FASTA, and JSON.

Accepts GenBank (.gb/.gbk) or SnapGene (.dna). SnapGene inputs are also written
as a normalized GenBank file. Features shorter than --min-length are skipped.
The JSON sidecar keeps product, note, and codon_start for hit annotation;
those strings do not fit in a FASTA header.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

# Feature types transferred onto integration-site hit windows.
FEATURE_TYPES = frozenset(
    {
        "CDS",
        "promoter",
        "terminator",
        "3'UTR",
        "misc_feature",
        "protein_bind",
    }
)


def detect_format(path: Path) -> str:
    ext = path.suffix.lower()
    if ext == ".dna":
        return "snapgene"
    if ext in {".gb", ".gbk", ".genbank"}:
        return "genbank"
    raise ValueError(f"Unsupported plasmid format for {path.name}: expected .gb/.gbk/.dna")


def sanitize_label(text: str) -> str:
    cleaned = re.sub(r"[^\w.\-]+", "_", text.strip())
    return cleaned.strip("_") or "feature"


def collapsed(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def first_qualifier(feat: SeqFeature, key: str) -> str | None:
    values = feat.qualifiers.get(key) or []
    if not values:
        return None
    text = collapsed(str(values[0]))
    return text or None


def joined_qualifier(feat: SeqFeature, key: str) -> str | None:
    values = feat.qualifiers.get(key) or []
    text = collapsed(" ".join(collapsed(str(value)) for value in values if str(value).strip()))
    return text or None


def raw_feature_label(feat: SeqFeature) -> str:
    for key in ("label", "gene", "product", "note", "name"):
        value = first_qualifier(feat, key)
        if value:
            return value
    return feat.type


def feature_label(feat: SeqFeature) -> str:
    return sanitize_label(raw_feature_label(feat))


def feature_metadata(feat: SeqFeature) -> dict[str, str]:
    """Qualifier text the hit annotator copies onto each transferred placement."""
    info = {"type": feat.type, "label": raw_feature_label(feat)}
    for key in ("product", "note"):
        value = joined_qualifier(feat, key)
        if value:
            info[key] = value
    codon_start = first_qualifier(feat, "codon_start")
    if codon_start:
        info["codon_start"] = codon_start
    return info


def feature_length(feat: SeqFeature) -> int:
    try:
        return len(feat.location)
    except Exception:
        return 0


def is_simple_location(feat: SeqFeature) -> bool:
    """Reject compound/join locations that cannot be extracted as one interval."""
    loc = feat.location
    if loc is None:
        return False
    if not isinstance(loc, FeatureLocation):
        # CompoundLocation and others
        parts = getattr(loc, "parts", None)
        return bool(parts) and len(parts) == 1 and isinstance(parts[0], FeatureLocation)
    return True


def simple_coords(feat: SeqFeature) -> tuple[int, int, int]:
    """Return 0-based start, end, and strand (+1/-1/0) for a simple feature."""
    loc = feat.location
    if not isinstance(loc, FeatureLocation):
        loc = loc.parts[0]
    start = int(loc.start)
    end = int(loc.end)
    strand = int(loc.strand) if loc.strand in (-1, 1) else 0
    return start, end, strand


def escape_gff_attr(value: str) -> str:
    return (
        value.replace(";", "%3B")
        .replace("=", "%3D")
        .replace("&", "%26")
        .replace(",", "%2C")
    )


def write_gff3(path: Path, seqid: str, features: list[tuple[str, SeqFeature]]) -> None:
    with path.open("w", encoding="utf-8") as out:
        out.write("##gff-version 3\n")
        for feat_id, feat in features:
            start, end, strand = simple_coords(feat)
            strand_char = "+" if strand >= 0 else "-"
            label = feature_label(feat)
            attrs = [
                f"ID={escape_gff_attr(feat_id)}",
                f"Name={escape_gff_attr(label)}",
            ]
            for key in ("product", "note", "gene"):
                if key in feat.qualifiers and feat.qualifiers[key]:
                    attrs.append(
                        f"{key}={escape_gff_attr(str(feat.qualifiers[key][0]))}"
                    )
            # GFF is 1-based inclusive
            out.write(
                "\t".join(
                    [
                        seqid,
                        "plasmid_extract_features",
                        feat.type,
                        str(start + 1),
                        str(end),
                        ".",
                        strand_char,
                        ".",
                        ";".join(attrs),
                    ]
                )
                + "\n"
            )


def extract_feature_seq(record: SeqRecord, feat: SeqFeature) -> str:
    try:
        return str(feat.extract(record.seq))
    except Exception as exc:
        raise ValueError(f"could not extract sequence: {exc}") from exc


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plasmid", type=Path, help="Input plasmid (.gb/.gbk/.dna)")
    parser.add_argument("--prefix", required=True, help="Output file prefix")
    parser.add_argument(
        "--min-length",
        type=int,
        default=8,
        help="Minimum feature length in bp (default: 8)",
    )
    args = parser.parse_args()

    in_format = detect_format(args.plasmid)
    record = SeqIO.read(str(args.plasmid), in_format)

    # Prefer a filesystem-friendly id for outputs / GFF seqid.
    seqid = sanitize_label(args.prefix) or sanitize_label(args.plasmid.stem)
    record.id = seqid
    record.name = seqid
    if not record.description:
        record.description = args.plasmid.stem

    gb_out = Path(f"{args.prefix}.gb")
    gff_out = Path(f"{args.prefix}.gff3")
    fa_out = Path(f"{args.prefix}.features.fasta")
    meta_out = Path(f"{args.prefix}.features.json")

    # Always write a GenBank copy (normalizes SnapGene input).
    SeqIO.write(record, str(gb_out), "genbank")

    kept: list[tuple[str, SeqFeature]] = []
    fasta_records: list[SeqRecord] = []
    metadata: dict[str, dict[str, str]] = {}
    skipped = 0

    for idx, feat in enumerate(record.features, start=1):
        if feat.type not in FEATURE_TYPES:
            continue
        if not is_simple_location(feat):
            print(
                f"WARN: skipping compound/unsupported location for {feat.type}",
                file=sys.stderr,
            )
            skipped += 1
            continue
        length = feature_length(feat)
        if length < args.min_length:
            print(
                f"WARN: skipping short {feat.type} ({length} bp < {args.min_length})",
                file=sys.stderr,
            )
            skipped += 1
            continue
        try:
            seq = extract_feature_seq(record, feat)
        except ValueError as exc:
            print(f"WARN: skipping {feat.type}: {exc}", file=sys.stderr)
            skipped += 1
            continue
        if len(seq) < args.min_length:
            skipped += 1
            continue

        feat_id = f"feat_{idx:04d}"
        label = feature_label(feat)
        kept.append((feat_id, feat))
        metadata[feat_id] = feature_metadata(feat)
        fasta_records.append(
            SeqRecord(
                seq=feat.extract(record.seq),
                id=feat_id,
                description=f"type={feat.type} label={label}",
            )
        )

    write_gff3(gff_out, seqid, kept)
    SeqIO.write(fasta_records, str(fa_out), "fasta")
    meta_out.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    print(
        f"Extracted {len(kept)} features (>= {args.min_length} bp); skipped {skipped}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
