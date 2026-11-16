#!/usr/bin/env python3
"""Build the Sniffles reference: unmodified host contigs plus plasmid contigs.

Plasmid sequence names get a plasmid_ prefix so they cannot collide with a host
contig. Host names are copied as written. A GenBank plasmid uses the LOCUS
name. This is the coordinate system for the structural-variant alignment, not
the extracted integration-site window.
"""

import argparse
import gzip
import re
import sys
from pathlib import Path

# A GenBank record ends with a line that is only "//". SnapGene headers
# contain https://, so splitting on every occurrence of those two characters
# would cut the record before ORIGIN.
GENBANK_RECORD_END = re.compile(r"(?m)^//\s*$")


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_fasta(text):
    name = None
    chunks = []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(chunks)
            # Keep the first token. The rest of the header is description.
            name = line[1:].split()[0]
            chunks = []
        elif name is not None:
            chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def parse_genbank(text):
    for record in GENBANK_RECORD_END.split(text):
        if "LOCUS" not in record or "ORIGIN" not in record:
            continue
        locus = None
        for line in record.splitlines():
            if line.startswith("LOCUS"):
                parts = line.split()
                if len(parts) < 2:
                    raise SystemExit(f"LOCUS line has no name: {line}")
                locus = parts[1]
                break
        if not locus:
            raise SystemExit("GenBank record has no LOCUS name")
        origin = record.split("ORIGIN", 1)[1]
        sequence = "".join(base for base in origin if base.isalpha())
        if not sequence:
            raise SystemExit(f"GenBank record {locus} has no sequence")
        yield locus, sequence


def read_sequences(path):
    with open_text(path) as handle:
        text = handle.read()
    stripped = text.lstrip()
    if not stripped:
        raise SystemExit(f"sequence file is empty: {path}")
    suffix = Path(path).name.lower()
    if suffix.endswith(".gz"):
        suffix = suffix[:-3]
    if stripped.startswith(">") or suffix.endswith((".fa", ".fasta", ".fna")):
        records = list(parse_fasta(text))
    elif stripped.startswith("LOCUS") or suffix.endswith((".gb", ".gbk", ".genbank")):
        records = list(parse_genbank(text))
    else:
        raise SystemExit(f"expected FASTA or GenBank: {path}")
    if not records:
        raise SystemExit(f"no sequences in {path}")
    return records


def plasmid_name(name):
    if name.startswith("plasmid_"):
        return name
    return f"plasmid_{name}"


def write_fasta(handle, records):
    for name, sequence in records:
        handle.write(f">{name}\n")
        for start in range(0, len(sequence), 60):
            handle.write(sequence[start : start + 60] + "\n")


def prepare(host_path, plasmid_path, output_path):
    host = read_sequences(host_path)
    plasmid = read_sequences(plasmid_path) if plasmid_path else []
    prefixed = [(plasmid_name(name), sequence) for name, sequence in plasmid]
    with open(output_path, "w") as handle:
        write_fasta(handle, host + prefixed)
    return [name for name, _sequence in host], [name for name, _sequence in prefixed]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", required=True, help="Host reference FASTA")
    parser.add_argument("--plasmid", help="Plasmid FASTA or GenBank. Omitted when the sample has none.")
    parser.add_argument("--output", required=True, help="Composite FASTA to write")
    args = parser.parse_args()
    host_names, plasmid_names = prepare(args.host, args.plasmid, args.output)
    print(
        f"host_contigs={len(host_names)} plasmid_contigs={len(plasmid_names)}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
