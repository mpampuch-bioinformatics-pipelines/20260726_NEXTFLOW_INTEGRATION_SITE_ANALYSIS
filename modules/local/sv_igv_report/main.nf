process SV_IGV_REPORT {
    tag "${meta.id}"
    label 'process_medium'

    // Same igv-reports image as the per-hit report. This HTML is a different
    // coordinate system: the host reference plus plasmid_ contigs.
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'docker://mpampuch/igv-reports_1.17.0_linux-amd64:d250c5633ca7d24f'
        : 'docker.io/mpampuch/igv-reports_1.17.0_linux-amd64:d250c5633ca7d24f'}"

    input:
    tuple val(meta), path(bam), path(bai), path(ref), path(vcf)

    output:
    tuple val(meta), path("${prefix}.igvreport.html"), emit: report
    tuple val("${task.process}"), val('igv-reports'), val('1.17.0'), topic: versions, emit: versions_igv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_structural_variants"
    """
    # Sniffles calls are on the composite reference, not the extracted hit.
    # A colon in a host accession is left as-is; only the per-hit report
    # rewrites :start-end suffixes.
    python3 - <<'PY'
import gzip
from pathlib import Path
src = Path("${vcf}")
text = gzip.open(src, "rt").read() if src.name.endswith(".gz") else src.read_text()
Path("sites.vcf").write_text(text)
PY
    records=\$(grep -c -v '^#' sites.vcf || true)
    if [[ "\$records" -eq 0 ]]; then
        printf '%s\\n' '<html><body><p>Sniffles2 wrote no structural-variant records for this sample. The VCF header is published next to this file.</p></body></html>' \\
            > ${prefix}.igvreport.html
    else
        samtools faidx "${ref}"
        create_report \\
            sites.vcf \\
            --fasta "${ref}" \\
            --tracks sites.vcf ${bam} \\
            --output ${prefix}.igvreport.html \\
            ${args}
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_structural_variants"
    """
    touch ${prefix}.igvreport.html
    """
}
