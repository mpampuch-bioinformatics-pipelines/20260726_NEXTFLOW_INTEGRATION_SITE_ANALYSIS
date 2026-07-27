process BLAST_MAKEBLASTDB {
    tag "${meta.id}_${meta.db_type}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_4'
        : 'quay.io/biocontainers/blast:2.16.0--h66d330f_4'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.n*"), emit: db
    tuple val("${task.process}"), val('blast'), eval("blastn -version | sed '1!d; s/blastn: //'"), topic: versions, emit: versions_blast

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // The DB basename encodes id + db_type so downstream blastn/blastdbcmd can
    // reconstruct it deterministically from the meta map.
    prefix = task.ext.prefix ?: "${meta.id}_${meta.db_type}"
    """
    makeblastdb \\
        -in ${fasta} \\
        -dbtype nucl \\
        -parse_seqids \\
        -out ${prefix} \\
        -title ${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${meta.db_type}"
    """
    touch ${prefix}.ndb ${prefix}.nhr ${prefix}.nin ${prefix}.not ${prefix}.nsq ${prefix}.ntf ${prefix}.nto
    """
}
