process DOWNSAMPLE_SEQKIT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
        : 'quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0'}"

    input:
    tuple val(meta), path(reads), path(genome_size)

    output:
    tuple val(meta), path("*.subsampled.fastq.gz")   , emit: reads
    tuple val(meta), path("*.downsample_summary.tsv"), emit: summary
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/seqkit v//'"), topic: versions, emit: versions_seqkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    def prefix          = task.ext.prefix ?: "${meta.id}"
    // Per-sample override falls back to the global default coverage target.
    def target_coverage = meta.target_coverage ?: params.target_coverage
    def seed            = params.seqkit_seed
    """
    # Genome size comes from a shared seqkit/stats TSV (sum_len), so the reference
    # is not scanned again for every sample that uses it.
    genome_size=\$(awk -F '\\t' 'NR==1 { for (i = 1; i <= NF; i++) if (\$i == "sum_len") col = i } NR==2 { print \$col }' ${genome_size})
    read_bases=\$(seqkit stats -T ${reads} | awk 'NR==2 {print \$5}')
    target_bases=\$(awk -v g="\$genome_size" -v c=${target_coverage} 'BEGIN { printf "%d", g * c }')
    fraction=\$(awk -v t="\$target_bases" -v b="\$read_bases" 'BEGIN { if (b <= 0) { print 0 } else { f = t / b; if (f > 1) f = 1; printf "%.6f", f } }')

    printf "sample\\tgenome_size\\tread_bases\\ttarget_coverage\\ttarget_bases\\tfraction\\n" > ${prefix}.downsample_summary.tsv
    printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" "${meta.id}" "\$genome_size" "\$read_bases" "${target_coverage}" "\$target_bases" "\$fraction" >> ${prefix}.downsample_summary.tsv

    seqkit sample \\
        ${args} \\
        -p "\$fraction" \\
        -s ${seed} \\
        ${reads} \\
        -o ${prefix}.subsampled.fastq.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.subsampled.fastq.gz
    touch ${prefix}.downsample_summary.tsv
    """
}
