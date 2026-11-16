#!/usr/bin/env python3
"""Re-annotate the curated hit windows and report missing plasmid features.

This is the check for 20261116_Annotation-NEEDS-FIXING. It is not part of the
default unit suite: it runs blastn against the real locus sequences.
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))

from annotation_compare import feature_keys, normalize_label

EXTRACT = ROOT / "bin" / "plasmid_extract_features.py"
ANNOTATE = ROOT / "bin" / "annotate_hit_from_blast.py"
EXAMPLES = ROOT / "20261116_Annotation-NEEDS-FIXING"
PLASMIDS = {
    "cmc103": ROOT / "PLASMIDS" / "CmC_103.gb",
    "cmc104": ROOT / "PLASMIDS" / "CmC_104.gb",
}


def run(cmd: list[str], **kwargs) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(cmd, check=False, capture_output=True, text=True, **kwargs)
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(cmd)}\n{result.stderr}"
        )
    return result


def plasmid_for(sample_name: str) -> Path:
    for prefix, path in PLASMIDS.items():
        if sample_name.startswith(prefix):
            return path
    raise ValueError(f"no plasmid mapping for {sample_name}")


def extract(plasmid: Path, dest: Path) -> Path:
    dest.mkdir(parents=True, exist_ok=True)
    prefix = dest / "plasmid"
    run(
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
    return prefix


def signatures(prefix: Path) -> set[tuple[str, str, int]]:
    meta = json.loads(prefix.with_suffix(".features.json").read_text(encoding="utf-8"))
    sigs: set[tuple[str, str, int]] = set()
    for rec in SeqIO.parse(str(prefix.with_suffix(".features.fasta")), "fasta"):
        info = meta[rec.id]
        sigs.add((info["type"], normalize_label(info["label"]), len(rec.seq)))
    return sigs


def annotate_window(prefix: Path, hit_gb: Path, dest: Path) -> Path:
    dest.mkdir(parents=True, exist_ok=True)
    record = SeqIO.read(str(hit_gb), "genbank")
    hit_fa = dest / "hit.fa"
    SeqIO.write(record, str(hit_fa), "fasta")
    db = dest / "hitdb"
    blast = dest / "feature_blast.txt"
    run(["makeblastdb", "-in", str(hit_fa), "-dbtype", "nucl", "-out", str(db)])
    run(
        [
            "blastn",
            "-task",
            "blastn-short",
            "-query",
            str(prefix.with_suffix(".features.fasta")),
            "-db",
            str(db),
            "-out",
            str(blast),
            "-outfmt",
            "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-evalue",
            "10",
            "-max_target_seqs",
            "500",
            "-num_threads",
            "4",
        ]
    )
    out = dest / "reannotated"
    run(
        [
            sys.executable,
            str(ANNOTATE),
            "--hit-fasta",
            str(hit_fa),
            "--features-fasta",
            str(prefix.with_suffix(".features.fasta")),
            "--features-meta",
            str(prefix.with_suffix(".features.json")),
            "--blast",
            str(blast),
            "--prefix",
            str(out),
            "--min-identity",
            "96",
        ]
    )
    return out.with_suffix(".gb")


def pair_loci(sample_dir: Path) -> tuple[list[tuple[Path, Path]], list[str]]:
    curated = {path.name: path for path in (sample_dir / "ground-truth").glob("*.gb")}
    pairs: list[tuple[Path, Path]] = []
    notes: list[str] = []
    used: set[str] = set()
    for predicted in sorted((sample_dir / "predicted").glob("*.gb")):
        match = next((name for name in curated if predicted.name.endswith(name)), None)
        if match is None:
            notes.append(f"unpaired predicted {predicted.name}")
            continue
        used.add(match)
        pairs.append((curated[match], predicted))
    for name in sorted(set(curated) - used):
        notes.append(f"unpaired curated {name}")
    return pairs, notes


def main() -> int:
    failures = 0
    with tempfile.TemporaryDirectory(prefix="plasmid-annot-") as tmp:
        work = Path(tmp)
        extracted: dict[Path, Path] = {}
        for sample_dir in sorted(p for p in EXAMPLES.iterdir() if p.is_dir()):
            plasmid = plasmid_for(sample_dir.name)
            if plasmid not in extracted:
                extracted[plasmid] = extract(plasmid, work / plasmid.stem)
            prefix = extracted[plasmid]
            sigs = signatures(prefix)
            pairs, notes = pair_loci(sample_dir)
            for note in notes:
                print(f"{sample_dir.name}: {note}")
            for curated, predicted in pairs:
                curated_rec = SeqIO.read(str(curated), "genbank")
                predicted_rec = SeqIO.read(str(predicted), "genbank")
                same_seq = str(curated_rec.seq) == str(predicted_rec.seq)
                dest = work / sample_dir.name / predicted.stem
                dest.mkdir(parents=True)
                # Curated and pipeline windows matched in the example set.
                # Coordinates below are on the curated record.
                annotated = annotate_window(prefix, curated, dest)
                curated_keys = feature_keys(curated)
                predicted_keys = feature_keys(annotated)
                expected = set()
                length_mismatch = []
                label_types = {(ftype, label) for ftype, label, _length in sigs}
                for key in curated_keys:
                    ftype, label, start, end = key
                    length = end - start + 1
                    if (ftype, label, length) in sigs:
                        expected.add(key)
                    elif (ftype, label) in label_types:
                        length_mismatch.append((key, length))
                missing = sorted(expected - predicted_keys)
                print(
                    f"{sample_dir.name}/{curated.name}: "
                    f"seq_match={same_seq} expected={len(expected)} "
                    f"missing={len(missing)} length_mismatch={len(length_mismatch)} "
                    f"annotated={len(predicted_keys)}"
                )
                for key in missing[:25]:
                    print(f"  missing {key[0]} {key[1]} {key[2]}..{key[3]}")
                if len(missing) > 25:
                    print(f"  ... {len(missing) - 25} more")
                for key, length in length_mismatch[:8]:
                    print(
                        f"  length differs from plasmid {key[0]} {key[1]} "
                        f"{key[2]}..{key[3]} ({length} bp)"
                    )
                if missing:
                    failures += len(missing)
                nos = ("terminator", "nos terminator", 50784, 51036)
                if curated.name.startswith("1_ptg000012l_441345-563981") and sample_dir.name.startswith(
                    "cmc103_33"
                ):
                    present = nos in predicted_keys
                    print(f"  cmc103_33 first NOS 50784..51036 present={present}")
                    if not present:
                        failures += 1
    if failures:
        print(f"FAILED: {failures} missing major placements")
        return 1
    print("OK: every plasmid-length curated feature was transferred")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
