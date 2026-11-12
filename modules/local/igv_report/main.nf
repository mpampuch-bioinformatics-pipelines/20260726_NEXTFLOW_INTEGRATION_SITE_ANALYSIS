process IGV_REPORT {
    tag "${meta.id}_${meta.db_type}_${meta.hit_id}"
    label 'process_medium'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://aangeloo/igv-reports-image:latest'
        : 'docker.io/aangeloo/igv-reports-image:latest'}"

    input:
    tuple val(meta), path(bam), path(bai), path(ref), path(gff), path(lifted, optional: true), path(loci)

    output:
    tuple val(meta), path("${prefix}.igvreport.html"), emit: report
    tuple val("${task.process}"), val('igv-reports'), val('1.0'), topic: versions, emit: versions_igv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def pad = params.igv_locus_pad
    prefix = task.ext.prefix ?: "${meta.id}_${meta.db_type}_${meta.hit_id}"

    """
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
    header=\$(grep "^>" reference.clean.fa | head -n 1 | cut -c2- | awk '{print \$1}')

    # One seqid, equal to the cleaned reference name. Any other value would
    # embed a track that never appears on the locus.
    gff_seqids=\$(awk '!/^#/ && \$1 != "" { print \$1 }' "${gff}" | sort -u)
    n_seqids=\$(printf '%s\\n' "\$gff_seqids" | awk 'NF { c++ } END { print c+0 }')
    if [ "\$n_seqids" -ne 1 ] || [ "\$gff_seqids" != "\$header" ]; then
        echo "GFF3 seqid '\$gff_seqids' does not match reference name '\$header'" >&2
        exit 1
    fi

    # Lifted genes keep the hit FASTA name. Rewrite the coordinate suffix the
    # same way as the FASTA header so the track sits on this locus. The
    # published lift file is not changed.
    tracks="${gff}"
    lifted_file="\${lifted:-}"
    if [[ -n "\$lifted_file" && -s "\$lifted_file" ]]; then
        case "\$lifted_file" in
            *.gtf) genes=reference_genes.gtf ;;
            *) genes=reference_genes.gff3 ;;
        esac
        awk 'BEGIN { FS = OFS = "\\t" }
            /^#/ { print; next }
            {
                if (\$1 ~ /:[0-9]+-[0-9]+\$/) {
                    seq = \$1
                    sub(/:[0-9]+-[0-9]+\$/, "", seq)
                    suffix = \$1
                    sub(/^.*:/, "", suffix)
                    \$1 = seq "_" suffix
                }
                print
            }
        ' "\$lifted_file" > "\$genes"
        gene_seqids=\$(awk '!/^#/ && \$1 != "" { print \$1 }' "\$genes" | sort -u)
        n_gene_seqids=\$(printf '%s\\n' "\$gene_seqids" | awk 'NF { c++ } END { print c+0 }')
        if [ "\$n_gene_seqids" -ne 1 ] || [ "\$gene_seqids" != "\$header" ]; then
            echo "Lifted annotation seqid '\$gene_seqids' does not match reference name '\$header'" >&2
            exit 1
        fi
        tracks="\$tracks \$genes"
    fi

    # The FASTA is the insertion plus --flank_size. Keep alignments that overlap
    # the insertion itself expanded by igv_locus_pad. A region query keeps a read
    # that spans the insertion or only crosses one edge; flank-only reads drop.
    interval=\$(igv_locus_interval.sh ${loci} "\$header" "\$len" ${pad})
    read -r keep_start keep_end _local_start _local_end <<< "\$interval"

    samtools view -b alignments.renamed.bam "\${header}:\${keep_start}-\${keep_end}" > locus.bam
    samtools index locus.bam

    # Cap relative to this window. A fraction of the unfiltered BAM would
    # discard almost every read that actually touches the insertion.
    count=\$(samtools view -c locus.bam)
    subsample=\$(echo \$count | awk '{ if (\$1 < 500) { print 1 } else { print 500 / \$1 } }')

    # BED is 0-based half-open and already includes the pad, so flanking
    # would open a margin that the filter just removed.
    bed_start=\$((keep_start - 1))
    echo -e "\${header}\t\${bed_start}\t\${keep_end}\tShowing \${subsample} fraction of \${count} alignments within ${pad} bp of the integration locus" > bedfile.bed

    create_report \\
        bedfile.bed \\
        --fasta reference.clean.fa \\
        --tracks \${tracks} locus.bam \\
        --output ${prefix}.igvreport.html \\
        --flanking 0 \\
        --subsample \${subsample} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${meta.db_type}_${meta.hit_id}"
    """
    touch ${prefix}.igvreport.html
    """
}
