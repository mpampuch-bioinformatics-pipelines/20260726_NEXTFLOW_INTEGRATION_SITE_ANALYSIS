process READ_INSPECTOR {
    tag "${meta.id}_${meta.db_type}"
    label 'process_high'

    // minimap2 2.30 and samtools/bgzip 1.23.1 from the installed nf-core
    // minimap2/align image. bgzip is on PATH in that image, so the script's
    // htslib step does not need a separate container.
    conda "${projectDir}/modules/nf-core/minimap2/align/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37671219cfd244eb9b33db9345d3543ffd83037419a1c57f4648aace493ec2c2/data'
        : 'community.wave.seqera.io/library/minimap2_samtools:b09096fc890429ce'}"

    input:
    tuple val(meta), path(reads), path(hits)

    output:
    tuple val(meta), path("read-inspector"), emit: inspector
    tuple val("${task.process}"), val('minimap2'), eval('minimap2 --version'), topic: versions, emit: versions_minimap2

    when:
    task.ext.when == null || task.ext.when

    script:
    // One directory for every hit of this sample and assembly. The session
    // files refer to reads.bam and contigs.fa.gz by relative path.
    """
    export PATH="/opt/conda/bin:\${PATH}"
    read_inspector.sh \\
        --outdir read-inspector \\
        --threads ${task.cpus} \\
        ${reads} \\
        ${hits}
    """

    stub:
    """
    mkdir -p read-inspector
    printf '>stub\\nACGT\\n' | gzip > read-inspector/contigs.fa.gz
    touch read-inspector/contigs.fa.gz.fai
    touch read-inspector/contigs.fa.gz.gzi
    touch read-inspector/reads.bam
    touch read-inspector/reads.bam.bai
    printf '%s\\n' '<?xml version="1.0"?><Session genome="contigs.fa.gz"><Resources><Resource path="reads.bam"/></Resources></Session>' > read-inspector/igv_session.xml
    printf '%s\\n' '{"assemblies":[{"sequence":{"adapter":{"fastaLocation":{"uri":"contigs.fa.gz"}}}}],"tracks":[{"adapter":{"bamLocation":{"uri":"reads.bam"}}}]}' > read-inspector/config.json
    """
}
