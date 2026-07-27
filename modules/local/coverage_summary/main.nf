process COVERAGE_SUMMARY {
    label 'process_low'

    // Custom image carries R + DT + sparkline needed by coverage_summary.R.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://aangeloo/nxf-tgs:latest'
        : 'docker.io/aangeloo/nxf-tgs:latest'}"

    input:
    path coverage_tsvs

    output:
    path "00-alignment-summary*", emit: summary
    tuple val("${task.process}"), val('coverage_summary'), val('1.0'), topic: versions, emit: versions_coverage_summary

    when:
    task.ext.when == null || task.ext.when

    script:
    // Prefer the Nextflow revision; fall back when running outside a git checkout.
    def commit = workflow.revision ?: workflow.scriptName ?: 'unknown'
    """
    coverage_summary.R "*.tsv" ${workflow.runName} ${commit}
    """

    stub:
    """
    touch 00-alignment-summary.html
    touch 00-alignment-summary.csv
    """
}
