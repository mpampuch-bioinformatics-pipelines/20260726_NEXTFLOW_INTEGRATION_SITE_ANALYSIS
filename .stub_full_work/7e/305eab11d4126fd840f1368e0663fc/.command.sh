#!/usr/bin/env bash -C -e -u -o pipefail
mkdir -p cmc-087-1-32_C.merolae_unpurged_raw_alternate
touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/report.tsv
touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/report.html
touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/report.pdf
touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/quast.log
touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/transposed_report.txt
touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/transposed_report.tex
touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/icarus.html
touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/report.tex
touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/report.txt

mkdir -p cmc-087-1-32_C.merolae_unpurged_raw_alternate/basic_stats
touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/basic_stats/cumulative_plot.pdf
touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/basic_stats/Nx_plot.pdf
touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/basic_stats/genome_GC_content_plot.pdf
touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/basic_stats/GC_content_plot.pdf

mkdir -p cmc-087-1-32_C.merolae_unpurged_raw_alternate/icarus_viewers
touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/icarus_viewers/contig_size_viewer.html

ln -s cmc-087-1-32_C.merolae_unpurged_raw_alternate/report.tsv cmc-087-1-32_C.merolae_unpurged_raw_alternate.tsv

if [ GCF_000091205.1_ASM9120v1_genomic.fna ]; then
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/basic_stats/NGx_plot.pdf
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/basic_stats/gc.icarus.txt

    mkdir -p cmc-087-1-32_C.merolae_unpurged_raw_alternate/aligned_stats
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/aligned_stats/NAx_plot.pdf
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/aligned_stats/NGAx_plot.pdf
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/aligned_stats/cumulative_plot.pdf

    mkdir -p cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/all_alignments_transcriptome.tsv
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/contigs_report_transcriptome.mis_contigs.info
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/contigs_report_transcriptome.stderr
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/contigs_report_transcriptome.stdout
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/contigs_report_transcriptome.unaligned.info
    mkdir -p cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/minimap_output
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/minimap_output/transcriptome.coords
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/minimap_output/transcriptome.coords.filtered
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/minimap_output/transcriptome.coords_tmp
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/minimap_output/transcriptome.sf
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/minimap_output/transcriptome.unaligned
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/minimap_output/transcriptome.used_snps
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/misassemblies_frcurve_plot.pdf
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/misassemblies_plot.pdf
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/misassemblies_report.tex
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/misassemblies_report.tsv
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/misassemblies_report.txt
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/transcriptome.mis_contigs.fa
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/transposed_report_misassemblies.tex
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/transposed_report_misassemblies.tsv
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/transposed_report_misassemblies.txt
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/unaligned_report.tex
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/unaligned_report.tsv
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/unaligned_report.txt

    mkdir -p cmc-087-1-32_C.merolae_unpurged_raw_alternate/genome_stats
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/genome_stats/genome_info.txt
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/genome_stats/transcriptome_gaps.txt
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/icarus_viewers/alignment_viewer.html

    ln -sf cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/misassemblies_report.tsv cmc-087-1-32_C.merolae_unpurged_raw_alternate_misassemblies.tsv
    ln -sf cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/unaligned_report.tsv cmc-087-1-32_C.merolae_unpurged_raw_alternate_unaligned.tsv
    ln -sf cmc-087-1-32_C.merolae_unpurged_raw_alternate/contigs_reports/all_alignments_transcriptome.tsv cmc-087-1-32_C.merolae_unpurged_raw_alternate_transcriptome.tsv

fi

if ([ GCF_000091205.1_ASM9120v1_genomic.fna ] && [  ]); then
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/genome_stats/features_cumulative_plot.pdf
    touch cmc-087-1-32_C.merolae_unpurged_raw_alternate/genome_stats/features_frcurve_plot.pdf
fi

# capture process environment
set +u
set +e
cd "$NXF_TASK_WORKDIR"

nxf_eval_cmd() {
    {
        IFS=$'\n' read -r -d '' "${1}";
        IFS=$'\n' read -r -d '' "${2}";
        (IFS=$'\n' read -r -d '' _ERRNO_; return ${_ERRNO_});
    } < <((printf '\0%s\0%d\0' "$(((({ shift 2; "${@}"; echo "${?}" 1>&3-; } | tr -d '\0' 1>&4-) 4>&2- 2>&1- | tr -d '\0' 1>&4-) 3>&1- | exit "$(cat)") 4>&1-)" "${?}" 1>&2) 2>&1)
}

echo '' > .command.env
#
nxf_eval_cmd STDOUT STDERR bash -c "quast.py --version 2>&1 | grep \"QUAST\" | sed 's/^.*QUAST v//; s/ .*$//'"
status=$?
if [ $status -eq 0 ]; then
  echo nxf_out_eval_12="$STDOUT" >> .command.env
  echo /nxf_out_eval_12/=exit:0 >> .command.env
else
  echo nxf_out_eval_12="$STDERR" >> .command.env
  echo /nxf_out_eval_12/=exit:$status >> .command.env
fi
