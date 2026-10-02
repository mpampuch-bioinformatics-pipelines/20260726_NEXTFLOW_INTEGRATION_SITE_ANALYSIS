#!/usr/bin/env bash -C -e -u -o pipefail
touch cmc_087.gb_6b0335e9_plasmid.gb
touch cmc_087.gb_6b0335e9_plasmid.gff3
printf ">feat_0001 type=CDS label=stub\nACGTACGTACGT\n" > cmc_087.gb_6b0335e9_plasmid.features.fasta
