process MODDOTPLOT {
    tag "${meta.id}_${meta.db_type}_${meta.hit_id}"
    label 'process_medium'

    // modDotPlot 1.0.0 (PyPI), frozen with Wave from containers/moddotplot.
    // Bioconda 0.9.9 cannot take --pairs. SingularityCE 3.9.7 cannot pull the
    // ORAS SIF. IBEX runs need Apptainer 1.4.5 on PATH (`mamba activate apptainer`).
    // oras:// bypasses NXF_SINGULARITY_CACHEDIR.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://docker.io/mpampuch/moddotplot_linux-amd64_sif:b20d8b6dee680bb8'
        : 'docker.io/mpampuch/moddotplot_linux-amd64:ede7fe8a21afa65b'}"

    input:
    tuple val(meta), path(pairs), path(reference), path(reference_fai), path(query), path(query_fai)

    output:
    tuple val(meta), path("moddotplot"), emit: moddotplot
    tuple val("${task.process}"), val('moddotplot'), eval("/opt/conda/bin/python -c 'import importlib.metadata as m; print(m.version(\"moddotplot\"))'"), topic: versions, emit: versions_moddotplot

    when:
    task.ext.when == null || task.ext.when

    script:
    // --pairs would otherwise compare every sequence in a multi-chromosome FASTA.
    // The index has to sit beside the FASTA; both are staged under their own names.
    """
    export PATH="/opt/conda/bin:\${PATH}"
    mapfile -t sequence_names < <(awk -F '\\t' '{ print \$1; print \$2 }' ${pairs} | sort -u)
    moddotplot static \\
        -f ${reference} ${query} \\
        -s "\${sequence_names[@]}" \\
        --pairs ${pairs} \\
        --compare-only \\
        -o moddotplot
    """

    stub:
    """
    mkdir -p moddotplot
    touch moddotplot/stub.txt
    """
}
