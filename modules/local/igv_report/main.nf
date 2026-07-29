process IGV_REPORT {
    tag "${meta.id}_${meta.db_type}_${meta.hit_id}"
    label 'process_medium'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://aangeloo/igv-reports-image:latest'
        : 'docker.io/aangeloo/igv-reports-image:latest'}"

    input:
    tuple val(meta), path(bam), path(bai), path(ref)

    output:
    tuple val(meta), path("${prefix}.igvreport.html"), emit: report
    tuple val("${task.process}"), val('igv-reports'), val('1.0'), topic: versions, emit: versions_igv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_${meta.db_type}_${meta.hit_id}"

    """
    # Subsample alignments above 500 so IGV HTML reports stay responsive.
    count=\$(samtools view -c ${bam})
    all_reads=\$(samtools flagstat ${bam} | grep 'primary\$' | cut -d" " -f1)
    subsample=\$(echo \$count | awk '{ if (\$1 < 500) { print 1 } else { print 500 / \$1 } }')

    # Create a FASTA with a cleaned sequence name (remove :start-end suffix)
    # Clean FASTA header:
    # 1. Remove trailing whitespace.
    # 2. Remove :start-end suffix from the sequence name.
    sed -E '
        /^>/ {
            s/[[:space:]]*\$//
            s/:[0-9]+-[0-9]+\$//
        }
    ' "${ref}" > reference.clean.fa

    samtools faidx reference.clean.fa

    len=\$(awk '{print \$2}' reference.clean.fa.fai | head -n 1)
    header=\$(grep "^>" reference.clean.fa | cut -c2- | head -n 1)

    echo -e "\${header}\t0\t\${len}\tShowing \${subsample} fraction of \${all_reads} alignments" > bedfile.bed

    create_report \\
        bedfile.bed \\
        --fasta reference.clean.fa \\
        --tracks ${bam} \\
        --output ${prefix}.igvreport.html \\
        --flanking 200 \\
        --subsample \${subsample} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${meta.db_type}_${meta.hit_id}"
    """
    touch ${prefix}.igvreport.html
    """
}