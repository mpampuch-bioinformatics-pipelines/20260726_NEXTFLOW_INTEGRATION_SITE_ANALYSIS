process JUNCTION_QC {
    tag "${meta.id}_${meta.db_type}"
    label 'process_low'

    // Python 3.12 and samtools 1.21, frozen with Wave from environment.yml.
    // The library tag is the freeze id, so the image does not float.
    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_samtools:f32393771b5fb0c6"

    input:
    tuple val(meta), path(inspector), path(loci)

    output:
    tuple val(meta), path("junction_qc.tsv"), emit: metrics
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed '1!d; s/samtools //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def anchor = params.junction_anchor_bp
    def window = params.junction_window_bp
    def pad = params.junction_depth_pad
    // PATH is set because the Wave image is micromamba and Nextflow's
    // singularity exec does not run the entrypoint.
    """
    export PATH="/opt/conda/bin:\${PATH}"
    junction_qc.py \\
        --sample '${meta.id}' \\
        --db-type '${meta.db_type}' \\
        --contigs ${inspector}/contigs.fa.gz \\
        --bam ${inspector}/reads.bam \\
        --loci ${loci} \\
        --anchor ${anchor} \\
        --window ${window} \\
        --depth-pad ${pad} \\
        --output junction_qc.tsv
    """

    stub:
    """
    printf '%s\\n' '# Flank-spanning reads, junction clips and splits, coverage, and reads per contig. This table does not choose trust or revisit.' > junction_qc.tsv
    printf '%s\\n' 'sample_id\tdb_type\tcontig\tlocus_id\tsubject\tlocal_start\tlocal_end\tseq_length\tleft_flank_bp\tright_flank_bp\tleft_flank_spanning\tright_flank_spanning\tleft_junction_clips\tright_junction_clips\tleft_junction_splits\tright_junction_splits\tleft_flank_mean_depth\tinsert_mean_depth\tright_flank_mean_depth\tleft_flank_zero_fraction\tinsert_zero_fraction\tright_flank_zero_fraction\tleft_junction_mean_depth\tright_junction_mean_depth\tn_reads\tnote' >> junction_qc.tsv
    printf '%s\\n' '${meta.id}\t${meta.db_type}\tstub\t1\tstub\t1\t1\t1\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t0\t' >> junction_qc.tsv
    """
}
