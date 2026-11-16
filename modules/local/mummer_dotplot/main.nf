process MUMMER_DOTPLOT {
    tag "${meta.id}_${meta.db_type}_${meta.hit_id}"
    label 'process_medium'

    // MUMmer 4.0.1 and gnuplot 5.4, frozen with Wave from containers/mummer4.
    // SingularityCE 3.9.7 cannot pull the ORAS SIF (Docker manifest media type).
    // IBEX runs need Apptainer 1.4.5 on PATH (`mamba activate apptainer`) so
    // `singularity` is Apptainer. An oras:// reference also bypasses
    // NXF_SINGULARITY_CACHEDIR.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://docker.io/mpampuch/mummer4_gnuplot_linux-amd64_sif:bc73a5fdeb8d3a78'
        : 'docker.io/mpampuch/mummer4_gnuplot_linux-amd64:01fef99d4c15605b'}"

    input:
    tuple val(meta), path(reference), path(query), path(best_fasta)

    output:
    tuple val(meta), path("${prefix}.dotplot.genome_mqc.png")          , emit: genome_png
    tuple val(meta), path("${prefix}.dotplot.chromosome.*_mqc.png")    , emit: chromosome_png
    tuple val("${task.process}"), val('mummer'), eval("/opt/conda/bin/nucmer --version 2>&1 | tail -1"), topic: versions, emit: versions_mummer

    when:
    task.ext.when == null || task.ext.when

    script:
    // mummerplot signals PID 1 after a non-x11 render (`kill 1, $child` with
    // $child left at 1). The task script is PID 1, so that SIGHUP is exit 129.
    // Hide gnuplot while mummerplot writes the script, then render it. gnuplot
    // 5.4 also rejects the "tiny" keyword on the png terminal line.
    prefix = task.ext.prefix ?: "${meta.id}_${meta.db_type}_${meta.hit_id}"
    """
    export PATH="/opt/conda/bin:\${PATH}"
    export GNUPLOT_DRIVER_DIR=/opt/conda/libexec/gnuplot/5.4
    render_mummerplot() {
        # This MUMmer 4 build runs gnuplot, then does `kill 1, \$child`.
        # For a png terminal \$child stays 1, so the signal hits PID 1.
        # Nextflow's task script is PID 1 and dies with SIGHUP (exit 129)
        # before it can fix the "png tiny" terminal. Hide gnuplot so
        # mummerplot only writes the script, then render that script.
        local bindir
        bindir=\$(mktemp -d)
        ln -s /opt/conda/bin/* "\$bindir"/
        rm -f "\$bindir/gnuplot"
        PATH="\$bindir:/usr/bin:/bin" mummerplot --fat --png --layout -p "\$2" "\$1"
        rm -rf "\$bindir"
        if [[ ! -s "\$2".gp ]]; then
            echo "mummerplot did not write \$2.gp" >&2
            exit 1
        fi
        sed -i -e 's/set terminal png tiny/set terminal png/' -e 's/set terminal qt.*/set terminal png/' "\$2".gp
        /opt/conda/bin/gnuplot "\$2".gp
    }

    nucmer --threads ${task.cpus} -p aln ${reference} ${query}
    if ! grep -q '^>' aln.delta; then
        echo "No nucmer alignment of the hit to the reference" >&2
        exit 1
    fi
    render_mummerplot aln.delta dotplot
    mv dotplot.png ${prefix}.dotplot.genome_mqc.png

    # The FASTA header is the accession chosen by PAF_HIT_REFERENCE. Put it
    # in the file name so the chromosome plot is not labelled "best_chromosome".
    chrom=\$(awk '/^>/ { sub(/^>/, ""); print \$1; exit }' ${best_fasta})
    safe=\$(printf '%s' "\$chrom" | tr -c 'A-Za-z0-9._-' '_')
    nucmer --threads ${task.cpus} -p aln.chromosome ${best_fasta} ${query}
    if ! grep -q '^>' aln.chromosome.delta; then
        echo "No nucmer alignment of the hit to the best chromosome" >&2
        exit 1
    fi
    render_mummerplot aln.chromosome.delta dotplot.chromosome
    mv dotplot.chromosome.png ${prefix}.dotplot.chromosome.\${safe}_mqc.png
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${meta.db_type}_${meta.hit_id}"
    """
    printf 'png' > ${prefix}.dotplot.genome_mqc.png
    printf 'png' > ${prefix}.dotplot.chromosome.chrom_long_mqc.png
    """
}
