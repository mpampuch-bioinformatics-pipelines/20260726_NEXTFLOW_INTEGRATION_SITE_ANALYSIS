process BLAST_BLASTN {
    tag "${meta.id}_${meta.db_type}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_4'
        : 'quay.io/biocontainers/blast:2.16.0--h66d330f_4'}"

    input:
    tuple val(meta), path(db), path(query)

    output:
    tuple val(meta), path("*.blast_results.txt"), emit: results
    tuple val("${task.process}"), val('blast'), eval("blastn -version | sed '1!d; s/blastn: //'"), topic: versions, emit: versions_blast

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.db_type}"
    // Matches the basename produced by BLAST_MAKEBLASTDB for this meta.
    def dbname = "${meta.id}_${meta.db_type}"
    def evalue = params.blast_evalue
    """
    blastn \\
        -query ${query} \\
        -db ${dbname} \\
        -out ${prefix}.blast_results.txt \\
        -evalue ${evalue} \\
        -num_threads ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.db_type}"
    """
    touch ${prefix}.blast_results.txt
    """
}
