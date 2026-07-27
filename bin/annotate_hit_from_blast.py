#!/usr/bin/env python3
"""Transfer plasmid features onto an extracted hit window using BLAST hits.

Reads blastn outfmt 6 (qseqid sseqid pident length mismatch gapopen qstart qend
sstart send evalue bitscore), keeps the best hit per feature at or above the
identity threshold, and writes GFF3 + GenBank annotations for the hit sequence.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


@dataclass
class BlastHit:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qlen: int | None = None

    @property
    def subject_start(self) -> int:
        return min(self.sstart, self.send)

    @property
    def subject_end(self) -> int:
        return max(self.sstart, self.send)

    @property
    def strand(self) -> int:
        return 1 if self.sstart <= self.send else -1

    @property
    def query_coverage(self) -> float:
        if not self.qlen or self.qlen <= 0:
            return 1.0
        return self.length / float(self.qlen)


def parse_fasta_headers(features_fasta: Path) -> dict[str, dict[str, str]]:
    """Map feat id -> {type, label} from description 'type=X label=Y'."""
    meta: dict[str, dict[str, str]] = {}
    for rec in SeqIO.parse(str(features_fasta), "fasta"):
        info = {"type": "misc_feature", "label": rec.id, "qlen": str(len(rec.seq))}
        for part in rec.description.split():
            if part.startswith("type="):
                info["type"] = part.split("=", 1)[1]
            elif part.startswith("label="):
                info["label"] = part.split("=", 1)[1]
        meta[rec.id] = info
    return meta


def parse_blast(path: Path, min_identity: float) -> list[BlastHit]:
    hits: list[BlastHit] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t") if "\t" in line else line.split()
            if len(cols) < 12:
                continue
            pident = float(cols[2])
            if pident < min_identity:
                continue
            hits.append(
                BlastHit(
                    qseqid=cols[0],
                    sseqid=cols[1],
                    pident=pident,
                    length=int(cols[3]),
                    qstart=int(cols[6]),
                    qend=int(cols[7]),
                    sstart=int(cols[8]),
                    send=int(cols[9]),
                    evalue=float(cols[10]),
                    bitscore=float(cols[11]),
                )
            )
    return hits


def best_hit_per_feature(
    hits: list[BlastHit],
    feature_meta: dict[str, dict[str, str]],
    min_query_coverage: float = 0.8,
) -> list[BlastHit]:
    by_query: dict[str, list[BlastHit]] = defaultdict(list)
    for hit in hits:
        meta = feature_meta.get(hit.qseqid)
        if meta and "qlen" in meta:
            hit.qlen = int(meta["qlen"])
        by_query[hit.qseqid].append(hit)

    selected: list[BlastHit] = []
    for _qid, group in by_query.items():
        # Prefer high bitscore, then identity, then coverage
        ranked = sorted(
            group,
            key=lambda h: (h.bitscore, h.pident, h.query_coverage),
            reverse=True,
        )
        best = ranked[0]
        if best.qlen and best.query_coverage < min_query_coverage:
            # Still accept if identity is high and length is short (blastn-short)
            if best.qlen >= 30:
                continue
        selected.append(best)
    return selected


def escape_gff_attr(value: str) -> str:
    return (
        value.replace(";", "%3B")
        .replace("=", "%3D")
        .replace("&", "%26")
        .replace(",", "%2C")
    )


def sanitize(text: str) -> str:
    return re.sub(r"[^\w.\-]+", "_", text.strip()).strip("_") or "feature"


def dedupe_hits(hits: list[BlastHit]) -> list[BlastHit]:
    """Drop identical transfers (same feature type/label/coords) from multi-HSP noise."""
    seen: set[tuple] = set()
    unique: list[BlastHit] = []
    for hit in hits:
        key = (hit.qseqid, hit.subject_start, hit.subject_end, hit.strand)
        if key in seen:
            continue
        seen.add(key)
        unique.append(hit)
    return unique


def write_gff3(
    path: Path,
    seqid: str,
    hits: list[BlastHit],
    feature_meta: dict[str, dict[str, str]],
) -> None:
    with path.open("w", encoding="utf-8") as out:
        out.write("##gff-version 3\n")
        for hit in hits:
            meta = feature_meta.get(hit.qseqid, {})
            ftype = meta.get("type", "misc_feature")
            label = meta.get("label", hit.qseqid)
            strand = "+" if hit.strand >= 0 else "-"
            attrs = [
                f"ID={escape_gff_attr(hit.qseqid)}",
                f"Name={escape_gff_attr(label)}",
                f"pident={hit.pident:.2f}",
                f"bitscore={hit.bitscore:.1f}",
            ]
            out.write(
                "\t".join(
                    [
                        seqid,
                        "plasmid_blast",
                        ftype,
                        str(hit.subject_start),
                        str(hit.subject_end),
                        f"{hit.pident:.2f}",
                        strand,
                        ".",
                        ";".join(attrs),
                    ]
                )
                + "\n"
            )


def write_genbank(
    path: Path,
    record: SeqRecord,
    hits: list[BlastHit],
    feature_meta: dict[str, dict[str, str]],
) -> None:
    features: list[SeqFeature] = [
        SeqFeature(
            FeatureLocation(0, len(record.seq), strand=1),
            type="source",
            qualifiers={"mol_type": ["other DNA"], "organism": ["integration hit"]},
        )
    ]
    for hit in hits:
        meta = feature_meta.get(hit.qseqid, {})
        ftype = meta.get("type", "misc_feature")
        label = meta.get("label", hit.qseqid)
        # Biopython FeatureLocation is 0-based, end exclusive; BLAST is 1-based inclusive
        start0 = hit.subject_start - 1
        end0 = hit.subject_end
        features.append(
            SeqFeature(
                FeatureLocation(start0, end0, strand=hit.strand),
                type=ftype,
                qualifiers={
                    "label": [label],
                    "note": [
                        f"transferred from plasmid feature {hit.qseqid}; "
                        f"pident={hit.pident:.2f}; bitscore={hit.bitscore:.1f}"
                    ],
                },
            )
        )
    out_rec = SeqRecord(
        seq=record.seq,
        id=sanitize(record.id),
        name=sanitize(record.id)[:16],
        description=record.description or "annotated integration hit",
        annotations={
            "molecule_type": "DNA",
            "topology": "linear",
            "data_file_division": "SYN",
        },
        features=features,
    )
    SeqIO.write(out_rec, str(path), "genbank")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hit-fasta", type=Path, required=True)
    parser.add_argument("--features-fasta", type=Path, required=True)
    parser.add_argument("--blast", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--min-identity", type=float, default=96.0)
    args = parser.parse_args()

    record = SeqIO.read(str(args.hit_fasta), "fasta")
    feature_meta = parse_fasta_headers(args.features_fasta)
    hits = parse_blast(args.blast, args.min_identity)
    selected = dedupe_hits(best_hit_per_feature(hits, feature_meta))

    seqid = sanitize(record.id.split()[0])
    gff_out = Path(f"{args.prefix}.gff3")
    gb_out = Path(f"{args.prefix}.gb")

    write_gff3(gff_out, seqid, selected, feature_meta)
    write_genbank(gb_out, record, selected, feature_meta)

    print(
        f"Annotated {len(selected)} features onto {seqid} "
        f"(from {len(hits)} BLAST rows >= {args.min_identity}% id)",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
