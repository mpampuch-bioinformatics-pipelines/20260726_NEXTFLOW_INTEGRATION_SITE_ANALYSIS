#!/usr/bin/env python3
"""Metrics for the four manual junction checks. No trust or revisit label."""

import tempfile
import unittest
from pathlib import Path

import junction_qc


CONTIG = "ptg000011l_1000-1199"
LENGTH = 200
# Insert is local 101-150. Left flank is 1-100. Right flank is 151-200.
LOCI = [
    "locus_id\tsubject\taln_start\taln_end\textract_start\textract_end\tn_members\tmembers",
    "1\tptg000011l\t1100\t1149\t1000\t1199\t1\tquery:1100-1149",
]


def sam(qname, flag, contig, pos, cigar, tags=""):
    query = cigar_query_length(cigar)
    optional = f"\t{tags}" if tags else ""
    return (
        f"{qname}\t{flag}\t{contig}\t{pos}\t60\t{cigar}\t*\t0\t0\t"
        f"{'A' * query}\t*{optional}"
    )


def cigar_query_length(cigar):
    total = 0
    for length, op in junction_qc.parse_cigar(cigar):
        if op in "MIS=X":
            total += length
    return total


def rows_for(sam_lines, depth_lines, contigs=None, loci=None, **kwargs):
    loci_rows = junction_qc.parse_loci(loci if loci is not None else LOCI)
    return junction_qc.build_rows(
        contigs if contigs is not None else [(CONTIG, LENGTH)],
        loci_rows,
        sam_lines,
        junction_qc.index_depth(depth_lines),
        sample="sample",
        db_type="combined",
        anchor=kwargs.get("anchor", 50),
        window=kwargs.get("window", 20),
        depth_pad=kwargs.get("depth_pad", 100),
    )


def depth_block(start, end, value, contig=CONTIG):
    return [f"{contig}\t{pos}\t{value}" for pos in range(start, end + 1)]


class JunctionQcTest(unittest.TestCase):
    def test_flank_spanning_needs_anchor_bases_on_both_sides(self):
        # 50 aligned bases on the flank and 50 on the insert. The short read
        # has only 10 bases past the junction.
        alignments = [
            sam("span_left", 0, CONTIG, 51, "100M"),
            sam("short_left", 0, CONTIG, 51, "60M"),
            sam("span_right", 0, CONTIG, 101, "100M"),
        ]
        row = rows_for(alignments, [])[0]
        self.assertEqual(row["left_flank_spanning"], "1")
        self.assertEqual(row["right_flank_spanning"], "1")
        self.assertEqual(row["local_start"], "101")
        self.assertEqual(row["local_end"], "150")
        self.assertNotIn("trust", row)
        self.assertNotIn("revisit", row)

    def test_clip_and_split_can_both_count_the_same_read(self):
        alignments = [
            sam("clip_only", 0, CONTIG, 101, "30S70M"),
            sam("both_read", 0, CONTIG, 101, "20S80M", f"SA:Z:{CONTIG},50,+,80M,60,0;"),
            # Supplementary piece starts on the right junction and has no clip.
            sam("split_right", 2048, CONTIG, 150, "50M", f"SA:Z:{CONTIG},50,+,80M,60,0;"),
        ]
        row = rows_for(alignments, [])[0]
        self.assertEqual(row["left_junction_clips"], "2")
        self.assertEqual(row["left_junction_splits"], "1")
        self.assertEqual(row["right_junction_clips"], "0")
        self.assertEqual(row["right_junction_splits"], "1")

    def test_insert_zero_run_does_not_lower_the_flank_mean(self):
        depth = []
        depth += depth_block(1, 100, 10)
        depth += depth_block(101, 125, 0)
        depth += depth_block(126, 150, 10)
        depth += depth_block(151, 200, 10)
        row = rows_for([], depth)[0]
        self.assertEqual(row["left_flank_mean_depth"], "10.0000")
        self.assertEqual(row["left_flank_zero_fraction"], "0.0000")
        self.assertEqual(row["insert_mean_depth"], "5.0000")
        self.assertEqual(row["insert_zero_fraction"], "0.5000")
        self.assertEqual(row["right_flank_mean_depth"], "10.0000")
        self.assertEqual(row["right_flank_zero_fraction"], "0.0000")

    def test_missing_left_flank_is_empty_not_zero(self):
        loci = [
            "locus_id\tsubject\taln_start\taln_end\textract_start\textract_end\tn_members\tmembers",
            "1\tptg000011l\t1000\t1049\t1000\t1199\t1\tquery:1000-1049",
        ]
        row = rows_for([], [], loci=loci)[0]
        self.assertEqual(row["local_start"], "1")
        self.assertEqual(row["left_flank_bp"], "")
        self.assertEqual(row["left_flank_spanning"], "")
        self.assertEqual(row["left_junction_clips"], "")
        self.assertEqual(row["left_junction_splits"], "")
        self.assertEqual(row["left_flank_mean_depth"], "")
        self.assertEqual(row["left_flank_zero_fraction"], "")
        self.assertEqual(row["left_junction_mean_depth"], "")
        self.assertEqual(row["right_flank_spanning"], "0")
        self.assertEqual(row["right_flank_mean_depth"], "0.0000")
        self.assertEqual(row["note"], "")

    def test_unmatched_contig_keeps_the_read_count(self):
        contigs = [(CONTIG, LENGTH), ("orphan", 10)]
        alignments = [sam("orphan_read", 0, "orphan", 1, "10M")]
        table = rows_for(alignments, [], contigs=contigs)
        self.assertEqual([row["contig"] for row in table], [CONTIG, "orphan"])
        orphan = table[1]
        self.assertEqual(orphan["n_reads"], "1")
        self.assertEqual(orphan["note"], "no_locus_row")
        self.assertEqual(orphan["locus_id"], "")
        self.assertEqual(orphan["left_flank_spanning"], "")
        self.assertEqual(orphan["insert_mean_depth"], "")
        self.assertEqual(orphan["seq_length"], "10")

    def test_report_has_no_decision_column(self):
        text = junction_qc.format_report(rows_for([], []))
        self.assertTrue(text.startswith("#"))
        self.assertIn("does not choose trust or revisit", text.splitlines()[0])
        header = next(line for line in text.splitlines() if not line.startswith("#"))
        fields = header.split("\t")
        self.assertEqual(fields, junction_qc.COLUMNS)
        for banned in ("trust", "revisit", "supported", "decision"):
            self.assertNotIn(banned, fields)

    def test_cli_reads_sam_and_depth_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / "contigs.fa"
            fasta.write_text(f">{CONTIG}\n{'A' * LENGTH}\n")
            loci = root / "loci.tsv"
            loci.write_text("\n".join(LOCI) + "\n")
            sam_path = root / "reads.sam"
            sam_path.write_text(sam("span_left", 0, CONTIG, 51, "100M") + "\n")
            depth_path = root / "depth.txt"
            depth_path.write_text("\n".join(depth_block(1, LENGTH, 1)) + "\n")
            out = root / "junction_qc.tsv"
            junction_qc.main(
                [
                    "--sample",
                    "sample",
                    "--db-type",
                    "combined",
                    "--contigs",
                    str(fasta),
                    "--loci",
                    str(loci),
                    "--sam",
                    str(sam_path),
                    "--depth",
                    str(depth_path),
                    "--output",
                    str(out),
                ]
            )
            text = out.read_text()
            self.assertIn("span_left", sam_path.read_text())
            data = next(line for line in text.splitlines() if line.startswith("sample\t"))
            fields = dict(zip(junction_qc.COLUMNS, data.split("\t")))
            self.assertEqual(fields["left_flank_spanning"], "1")
            self.assertEqual(fields["n_reads"], "1")


if __name__ == "__main__":
    unittest.main()
