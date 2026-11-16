process PLASMID_EXTRACT_FEATURES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/biopython:1.84'
        : 'quay.io/biocontainers/biopython:1.84'}"

    input:
    tuple val(meta), path(plasmid)

    output:
    tuple val(meta), path("${prefix}.gb")             , emit: gb
    tuple val(meta), path("${prefix}.gff3")           , emit: gff3
    tuple val(meta), path("${prefix}.features.fasta"), emit: features_fasta
    tuple val(meta), path("${prefix}.features.json"), emit: features_meta
    tuple val("${task.process}"), val('biopython'), val('1.84'), topic: versions, emit: versions_biopython

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}_plasmid"
    def minlen = params.plasmid_min_feature_length
    """
    plasmid_extract_features.py \\
        ${plasmid} \\
        --prefix ${prefix} \\
        --min-length ${minlen} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_plasmid"
    """
    touch ${prefix}.gb
    touch ${prefix}.gff3
    printf ">feat_0001 type=CDS label=stub\\nACGTACGTACGT\\n" > ${prefix}.features.fasta
    printf '%s\\n' '{"feat_0001":{"type":"CDS","label":"stub"}}' > ${prefix}.features.json
    """
}
