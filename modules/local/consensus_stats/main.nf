process CONSENSUS_STATS {
    tag "${meta.id}_${meta.db_type}_${meta.hit_id}"
    label 'process_low'

    // pegi3s image has Biopython plus basic Unix tools (awk) needed by run-average.awk
    conda "${moduleDir}/environment.yml"
    container 'docker.io/pegi3s/biopython:latest'

    input:
    tuple val(meta), path(consensus), path(coverage)

    output:
    tuple val(meta), path("${prefix}.coverage2.tsv"), emit: spark
    tuple val("${task.process}"), val('biopython'), val('1.84'), topic: versions, emit: versions_biopython

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_${meta.db_type}_${meta.hit_id}"
    """
    # Convert consensus FASTQ qualities to PHD, then subsample a rolling average for sparklines.
    fastq2phd.py -i ${consensus} -o ${prefix}.consensus.phd
    cat ${prefix}.consensus.phd | grep -v -e "BEGIN_*" -e "END*" > phd.file

    endpos=\$(cat phd.file | wc -l | awk '\$1=\$1' -)
    nth=\$(awk -v var="\$endpos" 'BEGIN { n = int(var / 300); if (n < 1) n = 1; print n }')
    run-average.awk col=2 size=10 phd.file | awk -v var="\$nth" 'NR % var == 0' - | tr '\\n' '|' > cons
    paste -d "\\t" ${coverage} cons > ${prefix}.coverage2.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${meta.db_type}_${meta.hit_id}"
    """
    cat ${coverage} > ${prefix}.coverage2.tsv
    echo -n "30|30|30|" >> ${prefix}.coverage2.tsv
    """
}
