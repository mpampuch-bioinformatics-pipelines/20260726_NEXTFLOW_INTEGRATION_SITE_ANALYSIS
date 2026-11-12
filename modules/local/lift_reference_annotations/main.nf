process LIFT_REFERENCE_ANNOTATIONS {
    tag "${meta.id}_${meta.db_type}_${meta.hit_id}"
    label 'process_medium'

    // Liftoff 1.6.1 and the minimap2 it depends on, frozen with Wave from
    // containers/liftoff. SingularityCE 3.9.7 cannot pull the ORAS SIF.
    // IBEX runs need Apptainer 1.4.5 on PATH (`mamba activate apptainer`).
    // oras:// bypasses NXF_SINGULARITY_CACHEDIR.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://docker.io/mpampuch/liftoff_linux-amd64_sif:da9cf7f71b08652c'
        : 'docker.io/mpampuch/liftoff_linux-amd64:bbc5403427a88eae'}"

    input:
    tuple val(meta), path(hit), path(reference), path(annotation)

    output:
    tuple val(meta), path("lifted.gff3"), optional: true, emit: lifted_gff3
    tuple val(meta), path("lifted.gtf") , optional: true, emit: lifted_gtf
    tuple val(meta), path("validation.tsv")         , emit: validation
    tuple val(meta), path("unmapped_features.txt")  , emit: unmapped
    tuple val("${task.process}"), val('liftoff'), eval("/opt/conda/bin/python -c 'import importlib.metadata as m; print(m.version(\"liftoff\"))'"), topic: versions, emit: versions_liftoff

    when:
    task.ext.when == null || task.ext.when

    script:
    // The bash script is the specification. A flank gene below coverage or
    // identity 0.5 makes it exit non-zero after writing validation.tsv.
    // A gene that spans the insertion is partial and does not fail the task.
    // No overlapping gene leaves validation.tsv and no lifted file.
    """
    export PATH="/opt/conda/bin:\${PATH}"
    lift_reference_annotations.sh \\
        --reference ${reference} \\
        --annotation ${annotation} \\
        --outdir lift_out \\
        --threads ${task.cpus} \\
        ${hit}
    name=\$(basename ${hit})
    name=\${name%.gz}
    name=\${name%.fa}
    name=\${name%.fasta}
    name=\${name%.fna}
    cp lift_out/\$name/validation.tsv validation.tsv
    cp lift_out/\$name/unmapped_features.txt unmapped_features.txt
    if [[ -s lift_out/\$name/lifted.gff3 ]]; then
        cp lift_out/\$name/lifted.gff3 lifted.gff3
    elif [[ -s lift_out/\$name/lifted.gtf ]]; then
        cp lift_out/\$name/lifted.gtf lifted.gtf
    fi
    """

    stub:
    """
    printf '##gff-version 3\\nchrom\\tstub\\tgene\\t1\\t3\\t.\\t+\\t.\\tID=near_gene\\n' > lifted.gff3
    printf 'gene_id\\tstatus\\tcoverage\\tsequence_id\\tvalid_orfs\\tnote\\nnear_gene\\tlifted\\t1\\t1\\t1\\t\\n' > validation.tsv
    touch unmapped_features.txt
    """
}
