process PREPARE_SV_REFERENCE {
    tag "${meta.id}"
    label 'process_single'

    // Same image as plasmid feature extraction. The script uses the Python
    // standard library, so a GenBank plasmid does not need a new container.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/biopython:1.84'
        : 'quay.io/biocontainers/biopython:1.84'}"

    input:
    tuple val(meta), path(host), path(plasmid)

    output:
    tuple val(meta), path("${prefix}.fasta"), emit: fasta
    tuple val("${task.process}"), val('biopython'), val('1.84'), topic: versions, emit: versions_biopython

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_sv_reference"
    def plasmid_arg = plasmid ? "--plasmid ${plasmid}" : ''
    """
    prepare_sv_reference.py \\
        --host ${host} \\
        --output ${prefix}.fasta \\
        ${plasmid_arg}
    """

    stub:
    // Stub-run stays on the host. The script is standard-library Python, so
    // the contig-prefix check does not need the container.
    prefix = task.ext.prefix ?: "${meta.id}_sv_reference"
    def plasmid_arg = plasmid ? "--plasmid ${plasmid}" : ''
    """
    prepare_sv_reference.py \\
        --host ${host} \\
        --output ${prefix}.fasta \\
        ${plasmid_arg}
    """
}
