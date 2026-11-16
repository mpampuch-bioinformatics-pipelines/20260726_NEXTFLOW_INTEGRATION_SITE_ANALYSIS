process PLASMID_ANNOTATE_HIT {
    tag "${meta.id}_${meta.db_type}_${meta.hit_id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/biopython:1.84'
        : 'quay.io/biocontainers/biopython:1.84'}"

    input:
    tuple val(meta), path(hit_fasta), path(features_fasta), path(blast_results), path(features_meta)

    output:
    tuple val(meta), path("${prefix}.gff3"), emit: gff3, optional: true
    tuple val(meta), path("${prefix}.gb")  , emit: gb  , optional: true
    tuple val("${task.process}"), val('biopython'), val('1.84'), topic: versions, emit: versions_biopython

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    prefix        = task.ext.prefix ?: "${meta.id}_${meta.db_type}_${meta.hit_id}"
    def min_ident = params.plasmid_min_identity
    """
    annotate_hit_from_blast.py \\
        --hit-fasta ${hit_fasta} \\
        --features-fasta ${features_fasta} \\
        --blast ${blast_results} \\
        --features-meta ${features_meta} \\
        --prefix ${prefix} \\
        --min-identity ${min_ident} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${meta.db_type}_${meta.hit_id}"
    """
    touch ${prefix}.gff3
    touch ${prefix}.gb
    """
}
