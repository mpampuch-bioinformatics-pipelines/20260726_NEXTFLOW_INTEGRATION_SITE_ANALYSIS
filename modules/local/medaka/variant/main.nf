process MEDAKA_VARIANT {
    tag "${meta.id}_${meta.db_type}_${meta.hit_id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/medaka:2.2.2--py312h3050eb1_0'
        : 'quay.io/biocontainers/medaka:2.2.2--py312h3050eb1_0'}"

    input:
    tuple val(meta), path(ref), path(reads)

    output:
    tuple val(meta), path("${prefix}.variants.vcf"), emit: vcf
    tuple val("${task.process}"), val('medaka'), val('1.11.3'), topic: versions, emit: versions_medaka

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}_${meta.db_type}_${meta.hit_id}"
    // Samplesheet may stage reads as a list; medaka_variant takes a single FASTQ.
    def input_reads = reads instanceof List ? reads[0] : reads
    """
    medaka_variant \\
        -i ${input_reads} \\
        -r ${ref} \\
        -t ${task.cpus} \\
        ${args}
    mv medaka/medaka.annotated.vcf ${prefix}.variants.vcf
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${meta.db_type}_${meta.hit_id}"
    """
    cat <<-END_VCF > ${prefix}.variants.vcf
    ##fileformat=VCFv4.1
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
    END_VCF
    """
}
