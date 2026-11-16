#!/usr/bin/env python3
"""Hit-window VCF contigs use the same colon rewrite as the IGV FASTA."""

import unittest

import igv_rewrite_vcf_contigs


VCF = """##fileformat=VCFv4.2
##contig=<ID=ptg000012l:441991-507307,length=65317>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
ptg000012l:441991-507307\t10\t.\tA\tC\t.\tPASS\tDP=8
"""


class RewriteVcfContigsTest(unittest.TestCase):
    def test_suffix_becomes_underscore_and_info_is_kept(self):
        rewritten = "".join(igv_rewrite_vcf_contigs.rewrite(VCF.splitlines(keepends=True)))
        self.assertIn("##contig=<ID=ptg000012l_441991-507307,length=65317>", rewritten)
        self.assertIn("ptg000012l_441991-507307\t10\t", rewritten)
        self.assertIn("DP=8", rewritten)
        self.assertNotIn("ptg000012l:441991-507307", rewritten)

    def test_name_without_a_coordinate_suffix_is_unchanged(self):
        self.assertEqual(igv_rewrite_vcf_contigs.rewrite_name("NC_000001.1"), "NC_000001.1")
        self.assertEqual(
            igv_rewrite_vcf_contigs.rewrite_name("ptg000012l:441991-507307"),
            "ptg000012l_441991-507307",
        )

    def test_expect_contig_rejects_a_different_sequence(self):
        with self.assertRaises(SystemExit):
            igv_rewrite_vcf_contigs.rewrite(
                VCF.splitlines(keepends=True),
                expect_contig="other_contig",
            )

    def test_empty_call_set_matches_the_cleaned_name(self):
        header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        rewritten = igv_rewrite_vcf_contigs.rewrite(
            header.splitlines(keepends=True),
            expect_contig="ptg000012l_441991-507307",
        )
        self.assertTrue(any(line.startswith("#CHROM") for line in rewritten))


if __name__ == "__main__":
    unittest.main()
