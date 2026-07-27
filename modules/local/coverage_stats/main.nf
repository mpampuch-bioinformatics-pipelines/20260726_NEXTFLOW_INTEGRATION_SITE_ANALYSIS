process COVERAGE_STATS {
    tag "${meta.id}_${meta.db_type}_${meta.hit_id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0'
        : 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${prefix}.coverage.tsv"), emit: coverage
    tuple val("${task.process}"), val('samtools'), val('1.21'), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_${meta.db_type}_${meta.hit_id}"
    // Sample label used as the first column of the coverage TSV (matches nxf-minimapper).
    def sample = prefix
    """
    # Subsample depth to ~300 points for sparkline plotting in the HTML summary.
    endpos=\$(samtools coverage -H ${bam} | cut -f3)
    nth=\$(awk -v var="\$endpos" 'BEGIN { n = int(var / 300); if (n < 1) n = 1; print n }')

    totalreads=\$(samtools flagstat ${bam} | grep 'primary\$' | cut -d" " -f1)
    mappedreads=\$(samtools flagstat ${bam} | grep 'primary mapped' | cut -d" " -f1)

    echo -e '${sample}'"\\t\$totalreads\\t\$mappedreads" > stats
    samtools coverage -H ${bam} > coverage
    samtools depth -aa ${bam} | awk -v var="\$nth" 'NR % var == 0' - | cut -f2,3 | tr '\\n' '|' | tr '\\t' ':' > depth
    paste -d "\\t" stats coverage depth > ${prefix}.coverage.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${meta.db_type}_${meta.hit_id}"
    """
    echo -e "${prefix}\\t100\\t80\\tref\\t1\\t100\\t80\\t100\\t100.0\\t10.0\\t30.0\\t60.0\\t1:10|50:20|" > ${prefix}.coverage.tsv
    """
}
