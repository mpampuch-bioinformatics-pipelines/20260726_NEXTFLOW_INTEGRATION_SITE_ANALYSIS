#!/usr/bin/env python3
"""Transfer plasmid features onto an extracted hit window using BLAST hits.

Reads blastn outfmt 6 (qseqid sseqid pident length mismatch gapopen qstart qend
sstart send evalue bitscore). Every placement at or above the identity and
coverage thresholds is kept, so tandem cassette copies are all annotated.
Overlapping HSPs of the same feature on the same strand are collapsed to the
stronger one. Writes GFF3 + GenBank for the hit sequence.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

# Oligos shorter than this often come back as a partial blastn-short HSP.
# They are still transferred when identity is high. Longer features must
# cover most of the query so a fragment is not stored as its own copy.
_SHORT_FEATURE_BP = 30


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
    placement_id: str = ""

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
    """Map feat id -> {type, label, qlen} from description 'type=X label=Y'."""
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


def load_feature_meta(
    features_fasta: Path,
    features_meta: Path | None,
) -> dict[str, dict[str, str]]:
    """Merge FASTA type/label/length with product, note, and codon_start from JSON."""
    meta = parse_fasta_headers(features_fasta)
    if features_meta is None:
        return meta
    payload = json.loads(features_meta.read_text(encoding="utf-8"))
    for feat_id, info in payload.items():
        if not isinstance(info, dict):
            continue
        slot = meta.setdefault(feat_id, {"type": "misc_feature", "label": feat_id})
        for key in ("type", "label", "product", "note", "codon_start"):
            value = info.get(key)
            if value:
                slot[key] = str(value)
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


def meets_query_coverage(hit: BlastHit, min_query_coverage: float) -> bool:
    if not hit.qlen or hit.qlen <= 0:
        return True
    return hit.query_coverage >= min_query_coverage


def _placement_rank(hit: BlastHit) -> tuple:
    return (hit.bitscore, hit.pident, hit.query_coverage, -hit.subject_start)


def same_feature_overlap(left: BlastHit, right: BlastHit) -> bool:
    """True when two HSPs of one feature describe the same locus on one strand."""
    if left.strand != right.strand:
        return False
    return left.subject_start <= right.subject_end and right.subject_start <= left.subject_end


def select_placements(
    hits: list[BlastHit],
    feature_meta: dict[str, dict[str, str]],
    min_query_coverage: float = 0.8,
) -> list[BlastHit]:
    """Keep every non-overlapping placement of each plasmid feature.

    A repeated cassette is several subject intervals of the same query. A
    weaker HSP that overlaps a better one is a fragment of that copy, not
    another copy. A short oligo keeps a single partial HSP only when no
    placement covers the query: with a loose BLAST e-value, 8 bp seeds of a
    12 bp linker would otherwise be stored as extra copies. Placement ids
    are numbered from the left of the hit.
    """
    by_query: dict[str, list[BlastHit]] = defaultdict(list)
    for hit in hits:
        meta = feature_meta.get(hit.qseqid)
        if meta and meta.get("qlen"):
            hit.qlen = int(meta["qlen"])
        by_query[hit.qseqid].append(hit)

    selected: list[BlastHit] = []
    for group in by_query.values():
        covered = [hit for hit in group if meets_query_coverage(hit, min_query_coverage)]
        if covered:
            pool = covered
        else:
            short = [
                hit
                for hit in group
                if hit.qlen is not None and 0 < hit.qlen < _SHORT_FEATURE_BP
            ]
            pool = [max(short, key=_placement_rank)] if short else []
        ranked = sorted(pool, key=_placement_rank, reverse=True)
        kept: list[BlastHit] = []
        for hit in ranked:
            if any(same_feature_overlap(hit, previous) for previous in kept):
                continue
            kept.append(hit)
        selected.extend(kept)

    selected.sort(key=lambda hit: (hit.subject_start, hit.subject_end, hit.qseqid, -hit.strand))
    counts: dict[str, int] = defaultdict(int)
    for hit in selected:
        counts[hit.qseqid] += 1
        hit.placement_id = f"{hit.qseqid}.{counts[hit.qseqid]}"
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
            placement = hit.placement_id or hit.qseqid
            attrs = [
                f"ID={escape_gff_attr(placement)}",
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
        qualifiers: dict[str, list[str]] = {"label": [label]}
        if meta.get("product"):
            qualifiers["product"] = [meta["product"]]
        if meta.get("codon_start"):
            qualifiers["codon_start"] = [meta["codon_start"]]
        notes: list[str] = []
        if meta.get("note"):
            notes.append(meta["note"])
        notes.append(
            f"transferred from plasmid feature {hit.qseqid}; "
            f"pident={hit.pident:.2f}; bitscore={hit.bitscore:.1f}"
        )
        qualifiers["note"] = notes
        features.append(
            SeqFeature(
                FeatureLocation(start0, end0, strand=hit.strand),
                type=ftype,
                qualifiers=qualifiers,
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
    parser.add_argument(
        "--features-meta",
        type=Path,
        default=None,
        help="JSON sidecar from plasmid_extract_features.py (product, note, codon_start)",
    )
    args = parser.parse_args()

    record = SeqIO.read(str(args.hit_fasta), "fasta")
    feature_meta = load_feature_meta(args.features_fasta, args.features_meta)
    hits = parse_blast(args.blast, args.min_identity)
    selected = select_placements(hits, feature_meta)

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
