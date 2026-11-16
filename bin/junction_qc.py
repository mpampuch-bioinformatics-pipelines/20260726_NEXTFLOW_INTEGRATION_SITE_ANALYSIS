#!/usr/bin/env python3
"""Junction QC metrics for manual review.

Counts flank-spanning reads, junction clips and splits, coverage across the
locus, and reads per candidate contig. These are the four manual checks.
The script does not choose trust or revisit.
"""

import argparse
import gzip
import re
import subprocess
import sys

# A read has to place this many aligned bases on the flank and on the insert.
DEFAULT_ANCHOR = 50
# A clip or split breakpoint counts when it falls within this many bases.
DEFAULT_WINDOW = 20
# Junction-local depth uses this many bases on each side of the boundary.
DEFAULT_DEPTH_PAD = 100

UNMAPPED = 0x4
SECONDARY = 0x100
SUPPLEMENTARY = 0x800

COLUMNS = [
    "sample_id",
    "db_type",
    "contig",
    "locus_id",
    "subject",
    "local_start",
    "local_end",
    "seq_length",
    "left_flank_bp",
    "right_flank_bp",
    "left_flank_spanning",
    "right_flank_spanning",
    "left_junction_clips",
    "right_junction_clips",
    "left_junction_splits",
    "right_junction_splits",
    "left_flank_mean_depth",
    "insert_mean_depth",
    "right_flank_mean_depth",
    "left_flank_zero_fraction",
    "insert_zero_fraction",
    "right_flank_zero_fraction",
    "left_junction_mean_depth",
    "right_junction_mean_depth",
    "n_reads",
    "note",
]

COMMENT = (
    "# Flank-spanning reads, junction clips and splits, coverage, and reads "
    "per contig. This table does not choose trust or revisit."
)

# Inspector rewrites the blastdbcmd colon (subject:start-end) to an underscore.
CONTIG_NAME = re.compile(r"^(.+)_([0-9]+)-([0-9]+)$")
CIGAR_OP = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_cigar(cigar):
    if cigar in ("", "*"):
        return []
    ops = [(int(length), op) for length, op in CIGAR_OP.findall(cigar)]
    consumed = sum(len(str(length)) + 1 for length, _op in ops)
    if consumed != len(cigar):
        raise ValueError(f"could not parse CIGAR {cigar}")
    return ops


def parse_loci(lines):
    rows = []
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if parts[0] == "locus_id":
            continue
        if len(parts) < 6:
            continue
        rows.append(
            {
                "locus_id": parts[0],
                "subject": parts[1],
                "aln_start": int(parts[2]),
                "aln_end": int(parts[3]),
                "extract_start": int(parts[4]),
                "extract_end": int(parts[5]),
            }
        )
    return rows


def contig_lengths(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    name = None
    length = 0
    order = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    order.append((name, length))
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line.strip())
    if name is not None:
        order.append((name, length))
    return order


def index_depth(lines):
    depth = {}
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        contig, pos, value = line.rstrip("\n").split("\t")[:3]
        depth.setdefault(contig, {})[int(pos)] = int(value)
    return depth


def match_locus(name, length, loci):
    """Return the insert span on this extracted sequence, or None.

    The same subject and extract_start rule as igv_locus_interval.sh. An exact
    extract_end wins. blastdbcmd can shorten the end when the requested range
    runs past the contig, so a longer table end still matches.
    """
    found = CONTIG_NAME.match(name)
    if not found:
        return None
    subject = found.group(1)
    header_start = int(found.group(2))
    header_end = int(found.group(3))
    if header_end < header_start:
        return None
    candidates = [
        row
        for row in loci
        if row["subject"] == subject
        and row["extract_start"] == header_start
        and row["extract_end"] >= header_end
    ]
    if not candidates:
        return None
    exact = [row for row in candidates if row["extract_end"] == header_end]
    if len(exact) > 1 or (not exact and len(candidates) > 1):
        raise ValueError(
            f"ambiguous locus rows for {subject} extract_start {header_start}"
        )
    chosen = exact[0] if exact else candidates[0]
    local_start = chosen["aln_start"] - chosen["extract_start"] + 1
    local_end = chosen["aln_end"] - chosen["extract_start"] + 1
    if local_start > local_end:
        local_start, local_end = local_end, local_start
    if local_end < 1 or local_start > length:
        return None
    # A truncated extraction still has a sequence. Clamp the span onto it.
    local_start = max(1, local_start)
    local_end = min(length, local_end)
    return chosen, local_start, local_end


def _aligned_sides(pos, cigar, junction, which):
    """Count M/=/X bases on the flank side and the insert side of a junction.

    Deletions are not aligned bases of the read. The left junction is the
    first insert base, so the flank is everything before it. The right
    junction is the last insert base, so the flank is everything after it.
    """
    flank = 0
    insert = 0
    ref = pos
    for length, op in parse_cigar(cigar):
        if op in "M=X":
            start = ref
            end = ref + length - 1
            if which == "left":
                flank_end = min(end, junction - 1)
                if start <= flank_end:
                    flank += flank_end - start + 1
                insert_start = max(start, junction)
                if insert_start <= end:
                    insert += end - insert_start + 1
            else:
                insert_end = min(end, junction)
                if start <= insert_end:
                    insert += insert_end - start + 1
                flank_start = max(start, junction + 1)
                if flank_start <= end:
                    flank += end - flank_start + 1
            ref += length
        elif op in "DN":
            ref += length
    return flank, insert


def alignment_end(pos, cigar):
    """Last reference base consumed by this alignment, or None if it has none."""
    ref = pos
    last = None
    for length, op in parse_cigar(cigar):
        if op in "MDN=X":
            last = ref + length - 1
            ref += length
    return last


def clip_breakpoints(pos, cigar):
    """Reference base where a soft or hard clip meets the alignment.

    SAM POS is already the first base after a leading clip, so a leading
    clip's breakpoint is POS. A trailing clip's breakpoint is the last
    aligned reference base.
    """
    ops = parse_cigar(cigar)
    if not ops:
        return []
    points = []
    leading = 0
    while leading < len(ops) and ops[leading][1] in "SH":
        leading += 1
    if leading and leading < len(ops):
        points.append(pos)
    trailing = len(ops) - 1
    while trailing >= 0 and ops[trailing][1] in "SH":
        trailing -= 1
    if trailing < len(ops) - 1 and trailing >= 0:
        end = alignment_end(pos, cigar)
        if end is not None:
            points.append(end)
    return points


def split_breakpoints(pos, cigar, flag, tags):
    """Where a split read breaks, if this record is one piece of a split.

    A clip alone is not a split. The record has to be supplementary or carry
    an SA tag. The breakpoint is that clip, and for a supplementary record
    also the start and end of its aligned piece: minimap2 often starts the
    second piece on the junction with no clip of its own.
    """
    supplementary = bool(flag & SUPPLEMENTARY)
    has_sa = any(tag.startswith("SA:Z:") for tag in tags)
    if not supplementary and not has_sa:
        return []
    points = clip_breakpoints(pos, cigar)
    if supplementary:
        points.append(pos)
        end = alignment_end(pos, cigar)
        if end is not None:
            points.append(end)
    return points


def _near(position, junction, window):
    return abs(position - junction) <= window


def _region_stats(depth_map, length, start, end):
    if start > end or start < 1:
        return None
    end = min(end, length)
    start = max(start, 1)
    if start > end:
        return None
    total = 0
    zeros = 0
    count = end - start + 1
    for pos in range(start, end + 1):
        value = depth_map.get(pos, 0)
        total += value
        if value == 0:
            zeros += 1
    return total / count, zeros / count


def _junction_mean(depth_map, length, last_base_before_cut, pad):
    """Mean depth in `pad` bases on each side of the cut.

    last_base_before_cut is the last base of the left-hand side. The cut for
    the left junction sits just before the insert. The cut for the right
    junction sits just after the insert.
    """
    if pad <= 0:
        return None
    start = last_base_before_cut - pad + 1
    end = last_base_before_cut + pad
    return _region_stats(depth_map, length, start, end)


def _fmt_float(stats, index):
    if stats is None:
        return ""
    return f"{stats[index]:.4f}"


def _fmt_mean(stats):
    if stats is None:
        return ""
    # _junction_mean returns the same (mean, zero_fraction) pair.
    return f"{stats[0]:.4f}"


def build_rows(contigs, loci, sam_lines, depth, sample, db_type, anchor, window, depth_pad):
    lengths = {name: length for name, length in contigs}
    placed = {}
    for name, length in contigs:
        try:
            matched = match_locus(name, length, loci)
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc
        placed[name] = matched

    primary = {name: set() for name, _length in contigs}
    spanning = {name: {"left": set(), "right": set()} for name, _length in contigs}
    clips = {name: {"left": set(), "right": set()} for name, _length in contigs}
    splits = {name: {"left": set(), "right": set()} for name, _length in contigs}

    for line in sam_lines:
        if not line.strip() or line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            continue
        qname = fields[0]
        flag = int(fields[1])
        rname = fields[2]
        if rname not in lengths or flag & UNMAPPED or flag & SECONDARY:
            continue
        pos = int(fields[3])
        cigar = fields[5]
        tags = fields[11:]
        is_primary = not (flag & SUPPLEMENTARY)
        if is_primary:
            primary[rname].add(qname)
        matched = placed[rname]
        if matched is None:
            continue
        _row, local_start, local_end = matched
        length = lengths[rname]
        has_left = local_start > 1
        has_right = local_end < length
        if is_primary and has_left:
            flank, insert = _aligned_sides(pos, cigar, local_start, "left")
            if flank >= anchor and insert >= anchor:
                spanning[rname]["left"].add(qname)
        if is_primary and has_right:
            flank, insert = _aligned_sides(pos, cigar, local_end, "right")
            if flank >= anchor and insert >= anchor:
                spanning[rname]["right"].add(qname)
        for point in clip_breakpoints(pos, cigar):
            if has_left and _near(point, local_start, window):
                clips[rname]["left"].add(qname)
            if has_right and _near(point, local_end, window):
                clips[rname]["right"].add(qname)
        for point in split_breakpoints(pos, cigar, flag, tags):
            if has_left and _near(point, local_start, window):
                splits[rname]["left"].add(qname)
            if has_right and _near(point, local_end, window):
                splits[rname]["right"].add(qname)

    rows = []
    for name, length in contigs:
        matched = placed[name]
        depth_map = depth.get(name, {})
        row = {column: "" for column in COLUMNS}
        row["sample_id"] = sample
        row["db_type"] = db_type
        row["contig"] = name
        row["seq_length"] = str(length)
        row["n_reads"] = str(len(primary[name]))
        if matched is None:
            row["note"] = "no_locus_row"
            rows.append(row)
            continue
        chosen, local_start, local_end = matched
        has_left = local_start > 1
        has_right = local_end < length
        row["locus_id"] = chosen["locus_id"]
        row["subject"] = chosen["subject"]
        row["local_start"] = str(local_start)
        row["local_end"] = str(local_end)
        row["note"] = ""
        if has_left:
            row["left_flank_bp"] = str(local_start - 1)
            row["left_flank_spanning"] = str(len(spanning[name]["left"]))
            row["left_junction_clips"] = str(len(clips[name]["left"]))
            row["left_junction_splits"] = str(len(splits[name]["left"]))
            left_stats = _region_stats(depth_map, length, 1, local_start - 1)
            row["left_flank_mean_depth"] = _fmt_float(left_stats, 0)
            row["left_flank_zero_fraction"] = _fmt_float(left_stats, 1)
            # The cut is the base before the first insert base.
            left_junction = _junction_mean(depth_map, length, local_start - 1, depth_pad)
            row["left_junction_mean_depth"] = _fmt_mean(left_junction)
        if has_right:
            row["right_flank_bp"] = str(length - local_end)
            row["right_flank_spanning"] = str(len(spanning[name]["right"]))
            row["right_junction_clips"] = str(len(clips[name]["right"]))
            row["right_junction_splits"] = str(len(splits[name]["right"]))
            right_stats = _region_stats(depth_map, length, local_end + 1, length)
            row["right_flank_mean_depth"] = _fmt_float(right_stats, 0)
            row["right_flank_zero_fraction"] = _fmt_float(right_stats, 1)
            right_junction = _junction_mean(depth_map, length, local_end, depth_pad)
            row["right_junction_mean_depth"] = _fmt_mean(right_junction)
        insert_stats = _region_stats(depth_map, length, local_start, local_end)
        row["insert_mean_depth"] = _fmt_float(insert_stats, 0)
        row["insert_zero_fraction"] = _fmt_float(insert_stats, 1)
        rows.append(row)
    return rows


def format_report(rows):
    lines = [COMMENT, "\t".join(COLUMNS)]
    for row in rows:
        lines.append("\t".join(row[column] for column in COLUMNS))
    return "\n".join(lines) + "\n"


def _read_text_lines(path):
    with open(path, "rt") as handle:
        return handle.readlines()


def _sam_and_depth_from_bam(bam):
    view = subprocess.Popen(
        ["samtools", "view", bam],
        stdout=subprocess.PIPE,
        text=True,
    )
    assert view.stdout is not None
    sam_lines = list(view.stdout)
    if view.wait() != 0:
        raise SystemExit(f"samtools view failed for {bam}")
    depth = subprocess.run(
        ["samtools", "depth", "-aa", bam],
        check=False,
        capture_output=True,
        text=True,
    )
    if depth.returncode != 0:
        raise SystemExit(depth.stderr.strip() or f"samtools depth failed for {bam}")
    return sam_lines, depth.stdout.splitlines()


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Write junction QC metrics. The table does not choose trust or revisit."
        )
    )
    parser.add_argument("--sample", required=True)
    parser.add_argument("--db-type", required=True)
    parser.add_argument("--contigs", required=True, help="FASTA or gzipped FASTA")
    parser.add_argument("--loci", required=True, help="loci.tsv or hit_intervals.tsv")
    parser.add_argument("--output", required=True)
    parser.add_argument("--bam", help="reads.bam from read inspector")
    parser.add_argument("--sam", help="samtools view text, used instead of --bam")
    parser.add_argument("--depth", help="samtools depth -aa text, used with --sam")
    parser.add_argument("--anchor", type=int, default=DEFAULT_ANCHOR)
    parser.add_argument("--window", type=int, default=DEFAULT_WINDOW)
    parser.add_argument("--depth-pad", type=int, default=DEFAULT_DEPTH_PAD)
    args = parser.parse_args(argv)

    if args.anchor < 1:
        raise SystemExit("--anchor must be at least 1")
    if args.window < 0 or args.depth_pad < 0:
        raise SystemExit("--window and --depth-pad must be zero or greater")

    if args.bam:
        sam_lines, depth_lines = _sam_and_depth_from_bam(args.bam)
    elif args.sam and args.depth:
        sam_lines = _read_text_lines(args.sam)
        depth_lines = _read_text_lines(args.depth)
    else:
        raise SystemExit("pass --bam, or both --sam and --depth")

    rows = build_rows(
        contig_lengths(args.contigs),
        parse_loci(_read_text_lines(args.loci)),
        sam_lines,
        index_depth(depth_lines),
        args.sample,
        args.db_type,
        args.anchor,
        args.window,
        args.depth_pad,
    )
    with open(args.output, "wt") as handle:
        handle.write(format_report(rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
