process BLAST_EXTRACT_HITS {
    tag "${meta.id}_${meta.db_type}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_4'
        : 'quay.io/biocontainers/blast:2.16.0--h66d330f_4'}"

    input:
    tuple val(meta), path(db), path(results)

    output:
    tuple val(meta), path("${prefix}.extracted.fa"), emit: fasta, optional: true
    tuple val(meta), path("all_hits/*.fa")         , emit: hits , optional: true
    tuple val("${task.process}"), val('blast'), eval("blastdbcmd -version | sed '1!d; s/blastdbcmd: //'"), topic: versions, emit: versions_blast

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix     = task.ext.prefix ?: "${meta.id}_${meta.db_type}"
    def dbname = "${meta.id}_${meta.db_type}"
    def flank  = params.flank_size
    """
    mkdir -p all_hits
    : > ${prefix}.extracted.fa
    count=0

    # For every BLAST hit (outfmt 6), extract the subject region padded by
    # +/- ${flank} bp of flanking sequence around the alignment coordinates.
    while read -r qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore; do
        [ -z "\$qseqid" ] && continue
        count=\$((count + 1))

        if [ "\$sstart" -lt "\$send" ]; then
            start=\$sstart
            end=\$send
        else
            start=\$send
            end=\$sstart
        fi

        ext_start=\$((start - ${flank}))
        [ "\$ext_start" -lt 1 ] && ext_start=1
        ext_end=\$((end + ${flank}))

        out="all_hits/\${count}_\${qseqid}_\${sseqid}_\${ext_start}-\${ext_end}.fa"
        blastdbcmd \\
            -db ${dbname} \\
            -entry "\$sseqid" \\
            -range "\${ext_start}-\${ext_end}" \\
            -outfmt "%f" \\
            -out "\$out" || { echo "WARN: extraction failed for \$sseqid" >&2; continue; }
        cat "\$out" >> ${prefix}.extracted.fa
    done < ${results}

    # Drop empty outputs so downstream mapping is only triggered when hits exist.
    [ -s ${prefix}.extracted.fa ] || rm -f ${prefix}.extracted.fa
    rmdir all_hits 2>/dev/null || true
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${meta.db_type}"
    """
    mkdir -p all_hits
    printf ">stub_hit\\nACGTACGTACGT\\n" > all_hits/1_stub.fa
    printf ">stub_hit\\nACGTACGTACGT\\n" > ${prefix}.extracted.fa
    """
}
