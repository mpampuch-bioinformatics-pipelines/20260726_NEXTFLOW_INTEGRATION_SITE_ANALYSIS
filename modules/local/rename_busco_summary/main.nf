process RENAME_BUSCO_SUMMARY {
    tag "${meta.id}_${meta.stage}_${meta.assembly}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:20.04'
        : 'nf-core/ubuntu:20.04'}"

    input:
    tuple val(meta), path(summaries)

    output:
    tuple val(meta), path("short_summary.*.${prefix}.txt"), emit: summary
    tuple val("${task.process}"), val('coreutils'), eval("cp --version | sed '1!d; s/.* //'"), topic: versions, emit: versions_coreutils

    when:
    task.ext.when == null || task.ext.when

    script:
    // Same id BUSCO_BUSCO uses, so the MultiQC row matches the assembly prefix.
    prefix = task.ext.prefix ?: "${meta.id}_${meta.stage}_${meta.assembly}"
    // Auto-lineage emits generic + specific short_summary files in one channel
    // item. Interpolating that collection as a single path concatenates names
    // and `cp` fails with "No such file". Copy each file on its own.
    """
    #region agent log
    names=\$(ls -1 short_summary*.txt 2>/dev/null | tr '\\n' ',')
    printf '{"sessionId":"b1da19","runId":"post-fix","hypothesisId":"A","location":"rename_busco_summary/main.nf:script","message":"staged short_summary files","data":{"names":"%s","prefix":"%s"},"timestamp":%s}\\n' "\$names" "${prefix}" "\$(date +%s)000" >> debug-b1da19.ndjson
    #endregion
    # Copy each staged summary. Do not interpolate the collection as one path.
    for summary in short_summary*.txt; do
        case "\$summary" in
            *".${prefix}.txt") continue ;;
        esac
        stem=\$(basename "\$summary" .txt)
        kind=\$(printf '%s\\n' "\$stem" | cut -d. -f2)
        lineage=\$(printf '%s\\n' "\$stem" | cut -d. -f3)
        cp -L "\$summary" "short_summary.\${kind}.\${lineage}.${prefix}.txt"
    done
    #region agent log
    out_names=\$(ls -1 short_summary.*.${prefix}.txt 2>/dev/null | tr '\\n' ',')
    printf '{"sessionId":"b1da19","runId":"post-fix","hypothesisId":"C","location":"rename_busco_summary/main.nf:after_cp","message":"renamed outputs","data":{"names":"%s"},"timestamp":%s}\\n' "\$out_names" "\$(date +%s)000" >> debug-b1da19.ndjson
    cat debug-b1da19.ndjson >> /Users/markpampuch/Downloads/tmp/20260726_NEXTFLOW_INTEGRATION_SITE_ANALYSIS/.cursor/debug-b1da19.log 2>/dev/null || true
    #endregion
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${meta.stage}_${meta.assembly}"
    """
    printf '# BUSCO version is: 0\\n' > short_summary.specific.lineage.${prefix}.txt
    """
}
