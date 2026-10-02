process BLAST_PREDICT_LOCI {
    tag "${meta.id}_${meta.db_type}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_4'
        : 'quay.io/biocontainers/blast:2.16.0--h66d330f_4'}"

    input:
    tuple val(meta), path(db), path(results)

    output:
    tuple val(meta), path("auto_loci_prediction/*.fa"), emit: loci, optional: true
    tuple val(meta), path("auto_loci_prediction/loci.tsv"), emit: summary
    tuple val("${task.process}"), val('blast'), eval("blastdbcmd -version | sed '1!d; s/blastdbcmd: //'"), topic: versions, emit: versions_blast

    when:
    task.ext.when == null || task.ext.when

    script:
    def dbname = "${meta.id}_${meta.db_type}"
    def flank  = params.flank_size
    def gap    = params.locus_merge_distance
    """
    mkdir -p auto_loci_prediction

    # One row per locus. Hits on the same contig within locus_merge_distance
    # bp are one insertion; the FASTA is the union of those alignments plus flanks.
    cluster_blast_loci.sh \\
        --merge-distance ${gap} \\
        --flank ${flank} \\
        ${results} \\
        > auto_loci_prediction/loci.tsv

    tail -n +2 auto_loci_prediction/loci.tsv | while IFS=\$'\\t' read -r locus_id subject _aln_start _aln_end extract_start extract_end _n_members _members; do
        [ -z "\$locus_id" ] && continue
        out="auto_loci_prediction/\${locus_id}_\${subject}_\${extract_start}-\${extract_end}.fa"
        blastdbcmd \\
            -db ${dbname} \\
            -entry "\$subject" \\
            -range "\${extract_start}-\${extract_end}" \\
            -outfmt "%f" \\
            -out "\$out" || { echo "WARN: locus extraction failed for \$subject" >&2; continue; }
    done

    # An empty glob is not a valid output when nothing clustered.
    find auto_loci_prediction -name '*.fa' -size 0 -delete
    """

    stub:
    """
    mkdir -p auto_loci_prediction
    printf ">stub_locus\\nACGTACGTACGT\\n" > auto_loci_prediction/1_stub_1-12.fa
    printf "locus_id\\tsubject\\taln_start\\taln_end\\textract_start\\textract_end\\tn_members\\tmembers\\n1\\tstub\\t1\\t12\\t1\\t12\\t1\\tstub:1-12\\n" > auto_loci_prediction/loci.tsv
    """
}
