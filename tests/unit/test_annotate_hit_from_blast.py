"""Each high-identity plasmid feature placement is transferred, not only the best one."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO

from annotate_hit_from_blast import BlastHit, select_placements

ROOT = Path(__file__).resolve().parents[2]
ANNOTATE = ROOT / "bin" / "annotate_hit_from_blast.py"


def _hit(
    *,
    qseqid: str = "feat_0001",
    pident: float = 100.0,
    length: int = 20,
    qstart: int = 1,
    qend: int = 20,
    sstart: int = 11,
    send: int = 30,
    bitscore: float = 40.0,
) -> BlastHit:
    return BlastHit(
        qseqid=qseqid,
        sseqid="hit",
        pident=pident,
        length=length,
        qstart=qstart,
        qend=qend,
        sstart=sstart,
        send=send,
        evalue=1e-8,
        bitscore=bitscore,
    )


def test_select_placements_keeps_disjoint_copies_and_drops_overlap() -> None:
    meta = {"feat_0001": {"qlen": "20", "type": "CDS", "label": "mVenus optimized"}}
    selected = select_placements(
        [
            _hit(sstart=11, send=30, length=20, bitscore=40),
            _hit(sstart=101, send=120, length=20, bitscore=39),
            _hit(sstart=11, send=25, length=15, qend=15, bitscore=10),
        ],
        meta,
    )
    assert [(hit.subject_start, hit.subject_end) for hit in selected] == [
        (11, 30),
        (101, 120),
    ]
    assert [hit.placement_id for hit in selected] == ["feat_0001.1", "feat_0001.2"]


def test_select_placements_drops_low_coverage_long_feature_keeps_short() -> None:
    long_meta = {"feat_0002": {"qlen": "100"}}
    assert (
        select_placements(
            [_hit(qseqid="feat_0002", length=40, qend=40, bitscore=80)],
            long_meta,
        )
        == []
    )
    short_meta = {"feat_0003": {"qlen": "12"}}
    short = select_placements(
        [_hit(qseqid="feat_0003", length=6, qend=6, sstart=5, send=10, bitscore=12)],
        short_meta,
    )
    assert len(short) == 1
    assert short[0].placement_id == "feat_0003.1"


def test_full_short_hits_suppress_distant_partial_seeds() -> None:
    meta = {"feat_0003": {"qlen": "12"}}
    selected = select_placements(
        [
            _hit(qseqid="feat_0003", length=12, qend=12, sstart=11, send=22, bitscore=24),
            _hit(qseqid="feat_0003", length=8, qend=8, sstart=50, send=57, bitscore=16),
        ],
        meta,
    )
    assert [(hit.subject_start, hit.subject_end) for hit in selected] == [(11, 22)]


def _write_hit(path: Path, length: int = 150) -> None:
    path.write_text(f">hit\n{'A' * length}\n", encoding="utf-8")


def _write_features(path: Path, seq: str, feat_id: str = "feat_0001") -> None:
    path.write_text(
        f">{feat_id} type=CDS label=mVenus_optimized\n{seq}\n",
        encoding="utf-8",
    )


def _blast_row(
    qseqid: str,
    pident: float,
    length: int,
    qend: int,
    sstart: int,
    send: int,
    bitscore: float,
) -> str:
    return "\t".join(
        [
            qseqid,
            "hit",
            f"{pident:.2f}",
            str(length),
            "0",
            "0",
            "1",
            str(qend),
            str(sstart),
            str(send),
            "1e-8",
            f"{bitscore:.1f}",
        ]
    )


def _run_annotate(
    tmp_path: Path,
    blast_rows: list[str],
    seq: str,
    *,
    meta: dict | None,
    min_identity: str = "96",
) -> Path:
    hit = tmp_path / "hit.fa"
    features = tmp_path / "features.fa"
    blast = tmp_path / "hits.txt"
    _write_hit(hit, max(150, len(seq) + 50))
    _write_features(features, seq)
    if meta is not None:
        (tmp_path / "features.json").write_text(json.dumps(meta), encoding="utf-8")
    blast.write_text("\n".join(blast_rows) + "\n", encoding="utf-8")
    prefix = tmp_path / "annotated"
    cmd = [
        sys.executable,
        str(ANNOTATE),
        "--hit-fasta",
        str(hit),
        "--features-fasta",
        str(features),
        "--blast",
        str(blast),
        "--prefix",
        str(prefix),
        "--min-identity",
        min_identity,
    ]
    if meta is not None:
        cmd.extend(["--features-meta", str(tmp_path / "features.json")])
    result = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if result.returncode != 0:
        raise AssertionError(f"{result.stderr}\n{result.stdout}")
    return prefix


def test_genbank_transfers_every_copy_and_qualifiers(tmp_path: Path) -> None:
    prefix = _run_annotate(
        tmp_path,
        [
            _blast_row("feat_0001", 100, 20, 20, 11, 30, 40),
            _blast_row("feat_0001", 100, 20, 20, 101, 120, 40),
            _blast_row("feat_0001", 100, 15, 15, 11, 25, 20),
        ],
        "A" * 20,
        meta={
            "feat_0001": {
                "type": "CDS",
                "label": "mVenus optimized",
                "product": "yellow fluorescent protein",
                "note": "codon optimized",
                "codon_start": "1",
            }
        },
    )
    record = SeqIO.read(str(prefix.with_suffix(".gb")), "genbank")
    cds = [feat for feat in record.features if feat.type == "CDS"]
    coords = [(int(feat.location.start) + 1, int(feat.location.end)) for feat in cds]
    assert coords == [(11, 30), (101, 120)]
    assert all(feat.qualifiers["label"][0] == "mVenus optimized" for feat in cds)
    assert all(feat.qualifiers["product"][0] == "yellow fluorescent protein" for feat in cds)
    assert all(feat.qualifiers["codon_start"][0] == "1" for feat in cds)
    assert all(
        any("codon optimized" in note for note in feat.qualifiers["note"]) for feat in cds
    )

    body = [
        line
        for line in prefix.with_suffix(".gff3").read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert [int(line.split("\t")[3]) for line in body] == [11, 101]
    attrs = [line.split("\t")[8] for line in body]
    assert "ID=feat_0001.1" in attrs[0]
    assert "ID=feat_0001.2" in attrs[1]


def test_cli_drops_below_identity_and_keeps_reverse_strand(tmp_path: Path) -> None:
    low = tmp_path / "low"
    low.mkdir()
    prefix = _run_annotate(
        low,
        [_blast_row("feat_0001", 90, 20, 20, 11, 30, 40)],
        "A" * 20,
        meta=None,
    )
    record = SeqIO.read(str(prefix.with_suffix(".gb")), "genbank")
    assert [feat.type for feat in record.features] == ["source"]

    reverse = tmp_path / "reverse"
    reverse.mkdir()
    prefix = _run_annotate(
        reverse,
        [_blast_row("feat_0001", 100, 20, 20, 30, 11, 40)],
        "A" * 20,
        meta=None,
    )
    record = SeqIO.read(str(prefix.with_suffix(".gb")), "genbank")
    cds = [feat for feat in record.features if feat.type == "CDS"]
    assert len(cds) == 1
    assert cds[0].location.strand == -1
    assert (int(cds[0].location.start) + 1, int(cds[0].location.end)) == (11, 30)
    assert cds[0].qualifiers["label"][0] == "mVenus_optimized"
