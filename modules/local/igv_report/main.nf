process IGV_REPORT {
    tag "${meta.id}_${meta.db_type}_${meta.hit_id}"
    label 'process_medium'

    // Custom image provides create_report (igv-reports) + samtools.
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://aangeloo/igv-reports-image:latest'
        : 'docker.io/aangeloo/igv-reports-image:latest'}"

    input:
    tuple val(meta), path(bam), path(bai), path(vcf), path(ref)

    output:
    tuple val(meta), path("${prefix}.igvreport.html"), emit: report
    tuple val("${task.process}"), val('igv-reports'), val('1.0'), topic: versions, emit: versions_igv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    prefix         = task.ext.prefix ?: "${meta.id}_${meta.db_type}_${meta.hit_id}"
    def vcf_track  = params.variants ? "${vcf}" : ""
    """
    # Subsample alignments above 500 so IGV HTML reports stay responsive.
    count=\$(samtools view -c ${bam})
    all_reads=\$(samtools flagstat ${bam} | grep 'primary\$' | cut -d" " -f1)
    subsample=\$(echo \$count | awk '{ if (\$1 < 500) { print 1 } else { print 500 / \$1 } }')

    samtools faidx ${ref}
    len=\$(awk '{print \$2}' ${ref}.fai | head -n 1)
    header=\$(grep ">" ${ref} | cut -c 2- | head -n 1)
    echo -e "\$header\\t0\\t\$len\\tShowing \$subsample fraction of \$all_reads alignments" > bedfile.bed

    if [ -n "${vcf_track}" ]; then
        vcf2bed.py -i ${vcf} -o vcf.bed
        head -n 10 vcf.bed >> bedfile.bed
    fi

    create_report \\
        bedfile.bed \\
        --fasta ${ref} \\
        --tracks ${vcf_track} ${bam} \\
        --output ${prefix}.igvreport.html \\
        --flanking 200 \\
        --subsample \$subsample \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${meta.db_type}_${meta.hit_id}"
    """
    touch ${prefix}.igvreport.html
    """
}
