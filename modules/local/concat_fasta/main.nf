process CONCAT_FASTA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:20.04'
        : 'nf-core/ubuntu:20.04'}"

    input:
    tuple val(meta), path(fastas, stageAs: "input*/*")

    output:
    tuple val(meta), path("${prefix}.fasta"), emit: fasta
    tuple val("${task.process}"), val('gzip'), eval("gzip --version | sed '1!d; s/.* //'"), topic: versions, emit: versions_gzip

    when:
    task.ext.when == null || task.ext.when

    script:
    // `db_type` (primary/alternate/combined) or a plain id keeps outputs unique.
    def suffix = meta.db_type ? "_${meta.db_type}" : ""
    prefix     = task.ext.prefix ?: "${meta.id}${suffix}"
    """
    # `zcat -f` transparently handles both gzipped and plain FASTA inputs.
    zcat -f ${fastas} > ${prefix}.fasta
    """

    stub:
    def suffix = meta.db_type ? "_${meta.db_type}" : ""
    prefix     = task.ext.prefix ?: "${meta.id}${suffix}"
    """
    touch ${prefix}.fasta
    """
}
