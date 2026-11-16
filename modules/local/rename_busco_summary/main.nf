process RENAME_BUSCO_SUMMARY {
    tag "${meta.id}_${meta.stage}_${meta.assembly}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:20.04'
        : 'nf-core/ubuntu:20.04'}"

    input:
    tuple val(meta), path(summary)

    output:
    tuple val(meta), path("short_summary.*.txt"), emit: summary
    tuple val("${task.process}"), val('coreutils'), eval("cp --version | sed '1!d; s/.* //'"), topic: versions, emit: versions_coreutils

    when:
    task.ext.when == null || task.ext.when

    script:
    // Same id BUSCO_BUSCO uses, so the MultiQC row matches the assembly prefix.
    prefix = task.ext.prefix ?: "${meta.id}_${meta.stage}_${meta.assembly}"
    """
    # BUSCO names every summary short_summary.specific.<lineage>.txt. MultiQC
    # only searches for short_summary* and keeps one file per basename, so a
    # primary and an alternate with the same lineage collapse to one sample.
    # Keep the short_summary prefix and put the haplotype id in front of .txt.
    stem=\$(basename "${summary}" .txt)
    kind=\$(printf '%s\\n' "\$stem" | cut -d. -f2)
    lineage=\$(printf '%s\\n' "\$stem" | cut -d. -f3)
    cp -L "${summary}" "short_summary.\${kind}.\${lineage}.${prefix}.txt"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${meta.stage}_${meta.assembly}"
    """
    printf '# BUSCO version is: 0\\n' > short_summary.specific.lineage.${prefix}.txt
    """
}
