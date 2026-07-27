process PLASMID_BLAST_FEATURES {
    tag "${meta.id}_${meta.db_type}_${meta.hit_id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_4'
        : 'quay.io/biocontainers/blast:2.16.0--h66d330f_4'}"

    input:
    tuple val(meta), path(hit_fasta), path(features_fasta)

    output:
    tuple val(meta), path(hit_fasta), path(features_fasta), path("${prefix}.feature_blast.txt"), emit: blast
    tuple val("${task.process}"), val('blast'), eval("blastn -version | sed '1!d; s/blastn: //'"), topic: versions, emit: versions_blast

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}_${meta.db_type}_${meta.hit_id}"
    def dbname = "${prefix}_hitdb"
    """
    if [ ! -s ${features_fasta} ]; then
        : > ${prefix}.feature_blast.txt
        exit 0
    fi

    makeblastdb \\
        -in ${hit_fasta} \\
        -dbtype nucl \\
        -out ${dbname}

    # blastn-short handles features down to ~8 bp.
    blastn \\
        -task blastn-short \\
        -query ${features_fasta} \\
        -db ${dbname} \\
        -out ${prefix}.feature_blast.txt \\
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \\
        -evalue ${params.blast_evalue} \\
        -num_threads ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${meta.db_type}_${meta.hit_id}"
    """
    touch ${prefix}.feature_blast.txt
    """
}
