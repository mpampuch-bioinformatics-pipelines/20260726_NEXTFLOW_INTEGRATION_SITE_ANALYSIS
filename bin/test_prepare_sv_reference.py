#!/usr/bin/env python3
"""Host contig names stay put. Plasmid contigs gain a plasmid_ prefix."""

import tempfile
import unittest
from pathlib import Path

import prepare_sv_reference


HOST = ">chrHost extra annotation\nACGTACGT\n"
PLASMID_FASTA = ">already\nTTTT\n>plasmid_kept\nAAAA\n"
GENBANK = """LOCUS       pMini                     8 bp    DNA     circular SYN 01-JAN-2026
DEFINITION  test plasmid
ACCESSION   .
JOURNAL   https://www.snapgene.com
FEATURES             Location/Qualifiers
ORIGIN
        1 acgtacgt
//
"""


class PrepareSvReferenceTest(unittest.TestCase):
    def test_host_names_unchanged_and_plasmid_prefixed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            host = root / "host.fna"
            plasmid = root / "plasmid.fa"
            output = root / "composite.fasta"
            host.write_text(HOST)
            plasmid.write_text(PLASMID_FASTA)
            host_names, plasmid_names = prepare_sv_reference.prepare(host, plasmid, output)
            self.assertEqual(host_names, ["chrHost"])
            self.assertEqual(plasmid_names, ["plasmid_already", "plasmid_kept"])
            text = output.read_text()
            self.assertIn(">chrHost\n", text)
            self.assertIn(">plasmid_already\n", text)
            self.assertIn(">plasmid_kept\n", text)
            self.assertNotIn(">already\n", text)

    def test_genbank_locus_becomes_plasmid_contig(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            host = root / "host.fna"
            plasmid = root / "plasmid.gb"
            output = root / "composite.fasta"
            host.write_text(HOST)
            plasmid.write_text(GENBANK)
            _host_names, plasmid_names = prepare_sv_reference.prepare(host, plasmid, output)
            self.assertEqual(plasmid_names, ["plasmid_pMini"])
            self.assertIn("ACGTACGT", output.read_text())

    def test_no_plasmid_keeps_only_the_host(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            host = root / "host.fna"
            output = root / "composite.fasta"
            host.write_text(HOST)
            host_names, plasmid_names = prepare_sv_reference.prepare(host, None, output)
            self.assertEqual(host_names, ["chrHost"])
            self.assertEqual(plasmid_names, [])
            self.assertEqual(output.read_text().count(">"), 1)


if __name__ == "__main__":
    unittest.main()
