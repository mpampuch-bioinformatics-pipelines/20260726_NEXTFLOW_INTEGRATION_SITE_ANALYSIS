process IGV_REPORT {
    tag "${meta.id}_${meta.db_type}_${meta.hit_id}"
    label 'process_medium'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://aangeloo/igv-reports-image:latest'
        : 'docker.io/aangeloo/igv-reports-image:latest'}"

    input:
    tuple val(meta), path(bam), path(bai), path(ref), path(gff)

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

    # A colon in the sequence name is parsed as an IGV locus. The GFF3 seqid
    # is the FASTA header with that colon replaced by an underscore
    # (ptg000012l:441991-507307 -> ptg000012l_441991-507307). Apply the same
    # replacement to the FASTA and the BAM @SQ SN so the annotation track,
    # alignments, and BED locus all name one sequence.
    sed -E '
        /^>/ {
            s/[[:space:]]*\$//
            s/:([0-9]+-[0-9]+)\$/_\\1/
        }
    ' "${ref}" > reference.clean.fa

    samtools faidx reference.clean.fa

    # Alignment records store a numeric ref id, so only @SQ SN needs to change.
    # [^:[:space:]] stops before the coordinate suffix; a greedy match would
    # swallow it and leave the original name unchanged.
    samtools view -H ${bam} \\
        | sed -E '/^@SQ/ s/(SN:[^:[:space:]]+):([0-9]+-[0-9]+)/\\1_\\2/' \\
        > header.sam
    samtools reheader header.sam ${bam} > alignments.renamed.bam
    samtools index alignments.renamed.bam

    len=\$(awk '{print \$2}' reference.clean.fa.fai | head -n 1)
    header=\$(grep "^>" reference.clean.fa | cut -c2- | head -n 1)

    # One seqid, equal to the cleaned reference name. Any other value would
    # embed a track that never appears on the locus.
    gff_seqids=\$(awk '!/^#/ && \$1 != "" { print \$1 }' "${gff}" | sort -u)
    n_seqids=\$(printf '%s\\n' "\$gff_seqids" | awk 'NF { c++ } END { print c+0 }')
    if [ "\$n_seqids" -ne 1 ] || [ "\$gff_seqids" != "\$header" ]; then
        echo "GFF3 seqid '\$gff_seqids' does not match reference name '\$header'" >&2
        exit 1
    fi

    echo -e "\${header}\t0\t\${len}\tShowing \${subsample} fraction of \${all_reads} alignments" > bedfile.bed

    create_report \\
        bedfile.bed \\
        --fasta reference.clean.fa \\
        --tracks "${gff}" alignments.renamed.bam \\
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
