#!/usr/bin/env python3
"""Rewrite hit-window contig names so IGV does not parse a colon as a locus.

The per-hit FASTA and BAM use the same substitution: a trailing
:start-end suffix becomes _start-end. Apply it to ##contig IDs and to the
CHROM column. INFO and sample fields are left alone.

With --expect-contig, every rewritten contig must be that one name. A header
with no contig and no records is an empty call set and is allowed.
"""

import argparse
import gzip
import re
import sys

# Same boundary as the FASTA sed in IGV_REPORT: the coordinate suffix is the
# last :digits-digits token, not an earlier colon in the sequence name.
CONTIG_SUFFIX = re.compile(r"^(?P<name>.+):(?P<bounds>[0-9]+-[0-9]+)$")
CONTIG_HEADER = re.compile(r"(ID=)([^,>\s]+)")


def rewrite_name(name):
    match = CONTIG_SUFFIX.match(name)
    if not match:
        return name
    return f"{match.group('name')}_{match.group('bounds')}"


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def rewrite(lines, expect_contig=None):
    contigs = set()
    output = []
    for line in lines:
        if line.startswith("##contig="):
            def replace(match):
                rewritten = rewrite_name(match.group(2))
                contigs.add(rewritten)
                return match.group(1) + rewritten

            output.append(CONTIG_HEADER.sub(replace, line, count=1))
            continue
        if line.startswith("#"):
            output.append(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if not fields or fields == [""]:
            continue
        fields[0] = rewrite_name(fields[0])
        contigs.add(fields[0])
        output.append("\t".join(fields) + "\n")
    if expect_contig is not None and contigs and contigs != {expect_contig}:
        found = ", ".join(sorted(contigs))
        raise SystemExit(
            f"VCF contig(s) {found!r} do not match reference name {expect_contig!r}"
        )
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_vcf")
    parser.add_argument("output_vcf")
    parser.add_argument("--expect-contig", help="Cleaned FASTA sequence name this VCF must use")
    args = parser.parse_args()
    with open_text(args.input_vcf) as handle:
        rewritten = rewrite(handle.readlines(), args.expect_contig)
    with open(args.output_vcf, "w") as handle:
        handle.writelines(rewritten)


if __name__ == "__main__":
    main()
