process MEDAKA_CONSENSUS {
    tag "${meta.id}_${meta.db_type}_${meta.hit_id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/medaka:1.11.3--py39h05d5c5f_0'
        : 'quay.io/biocontainers/medaka:1.11.3--py39h05d5c5f_0'}"

    input:
    tuple val(meta), path(ref), path(reads)

    output:
    tuple val(meta), path("${prefix}.consensus.fastq"), emit: consensus
    tuple val("${task.process}"), val('medaka'), val('1.11.3'), topic: versions, emit: versions_medaka

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}_${meta.db_type}_${meta.hit_id}"
    // Samplesheet may stage reads as a list; medaka_consensus takes a single FASTQ.
    def input_reads = reads instanceof List ? reads[0] : reads
    """
    # Fill uncovered gaps with N and emit FASTQ qualities for downstream sparkline stats.
    medaka_consensus \\
        -r N \\
        -q \\
        -i ${input_reads} \\
        -d ${ref} \\
        -t ${task.cpus} \\
        ${args}
    mv medaka/consensus.fastq ${prefix}.consensus.fastq
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${meta.db_type}_${meta.hit_id}"
    """
    printf "@stub\\nACGTACGT\\n+\\nIIIIIIII\\n" > ${prefix}.consensus.fastq
    """
}
