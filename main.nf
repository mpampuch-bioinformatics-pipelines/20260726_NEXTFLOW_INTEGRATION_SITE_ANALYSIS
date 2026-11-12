#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/integrationsiteanalysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Genomic integration-site analysis from long reads.
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INTEGRATION_SITE_ANALYSIS } from './workflows/integration_site_analysis'
include { PIPELINE_INITIALISATION   } from './subworkflows/local/utils_nfcore_pipeline_pipeline'
include { PIPELINE_COMPLETION       } from './subworkflows/local/utils_nfcore_pipeline_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline
//
workflow NFCORE_INTEGRATIONSITEANALYSIS {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    INTEGRATION_SITE_ANALYSIS (
        samplesheet
    )
    emit:
    multiqc_report      = INTEGRATION_SITE_ANALYSIS.out.multiqc_report // channel: /path/to/multiqc_report.html
    multiqc_data        = INTEGRATION_SITE_ANALYSIS.out.multiqc_data
    multiqc_plots       = INTEGRATION_SITE_ANALYSIS.out.multiqc_plots
    versions            = INTEGRATION_SITE_ANALYSIS.out.versions
    read_stats          = INTEGRATION_SITE_ANALYSIS.out.read_stats
    genomescope         = INTEGRATION_SITE_ANALYSIS.out.genomescope
    genomescope_plots   = INTEGRATION_SITE_ANALYSIS.out.genomescope_plots
    downsample_summary  = INTEGRATION_SITE_ANALYSIS.out.downsample_summary
    draft_assemblies    = INTEGRATION_SITE_ANALYSIS.out.draft_assemblies
    asmstats_raw        = INTEGRATION_SITE_ANALYSIS.out.asmstats_raw
    gfastats_raw        = INTEGRATION_SITE_ANALYSIS.out.gfastats_raw
    busco_raw           = INTEGRATION_SITE_ANALYSIS.out.busco_raw
    quast_raw           = INTEGRATION_SITE_ANALYSIS.out.quast_raw
    asmstats_purged     = INTEGRATION_SITE_ANALYSIS.out.asmstats_purged
    gfastats_purged     = INTEGRATION_SITE_ANALYSIS.out.gfastats_purged
    busco_purged        = INTEGRATION_SITE_ANALYSIS.out.busco_purged
    quast_purged        = INTEGRATION_SITE_ANALYSIS.out.quast_purged
    purged_assemblies   = INTEGRATION_SITE_ANALYSIS.out.purged_assemblies
    purgedups_files     = INTEGRATION_SITE_ANALYSIS.out.purgedups_files
    blast_results       = INTEGRATION_SITE_ANALYSIS.out.blast_results
    extracted_hits      = INTEGRATION_SITE_ANALYSIS.out.extracted_hits
    plasmid_gb          = INTEGRATION_SITE_ANALYSIS.out.plasmid_gb
    plasmid_gff3        = INTEGRATION_SITE_ANALYSIS.out.plasmid_gff3
    plasmid_features    = INTEGRATION_SITE_ANALYSIS.out.plasmid_features
    hit_annotations_gff = INTEGRATION_SITE_ANALYSIS.out.hit_annotations_gff
    hit_annotations_gb  = INTEGRATION_SITE_ANALYSIS.out.hit_annotations_gb
    mapped_bam          = INTEGRATION_SITE_ANALYSIS.out.mapped_bam
    mapped_bai          = INTEGRATION_SITE_ANALYSIS.out.mapped_bai
    coverage_summary    = INTEGRATION_SITE_ANALYSIS.out.coverage_summary
    igv_report          = INTEGRATION_SITE_ANALYSIS.out.igv_report
    dotplot_paf         = INTEGRATION_SITE_ANALYSIS.out.dotplot_paf
    dotplot_gaps        = INTEGRATION_SITE_ANALYSIS.out.dotplot_gaps
    dotplot_best        = INTEGRATION_SITE_ANALYSIS.out.dotplot_best
    dotplot_best_fasta  = INTEGRATION_SITE_ANALYSIS.out.dotplot_best_fasta
    dotplot_genome      = INTEGRATION_SITE_ANALYSIS.out.dotplot_genome
    dotplot_chromosome  = INTEGRATION_SITE_ANALYSIS.out.dotplot_chromosome
    dotplot_moddotplot  = INTEGRATION_SITE_ANALYSIS.out.dotplot_moddotplot
    lifted_annotation   = INTEGRATION_SITE_ANALYSIS.out.lifted_annotation
    lifted_validation   = INTEGRATION_SITE_ANALYSIS.out.lifted_validation
    lifted_unmapped     = INTEGRATION_SITE_ANALYSIS.out.lifted_unmapped
    read_inspector      = INTEGRATION_SITE_ANALYSIS.out.read_inspector
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_INTEGRATIONSITEANALYSIS (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_INTEGRATIONSITEANALYSIS.out.multiqc_report
    )

    publish:
    multiqc_report      = NFCORE_INTEGRATIONSITEANALYSIS.out.multiqc_report
    multiqc_data        = NFCORE_INTEGRATIONSITEANALYSIS.out.multiqc_data
    multiqc_plots       = NFCORE_INTEGRATIONSITEANALYSIS.out.multiqc_plots
    versions            = NFCORE_INTEGRATIONSITEANALYSIS.out.versions
    read_stats          = NFCORE_INTEGRATIONSITEANALYSIS.out.read_stats
    genomescope         = NFCORE_INTEGRATIONSITEANALYSIS.out.genomescope
    genomescope_plots   = NFCORE_INTEGRATIONSITEANALYSIS.out.genomescope_plots
    downsample_summary  = NFCORE_INTEGRATIONSITEANALYSIS.out.downsample_summary
    draft_assemblies    = NFCORE_INTEGRATIONSITEANALYSIS.out.draft_assemblies
    asmstats_raw        = NFCORE_INTEGRATIONSITEANALYSIS.out.asmstats_raw
    gfastats_raw        = NFCORE_INTEGRATIONSITEANALYSIS.out.gfastats_raw
    busco_raw           = NFCORE_INTEGRATIONSITEANALYSIS.out.busco_raw
    quast_raw           = NFCORE_INTEGRATIONSITEANALYSIS.out.quast_raw
    asmstats_purged     = NFCORE_INTEGRATIONSITEANALYSIS.out.asmstats_purged
    gfastats_purged     = NFCORE_INTEGRATIONSITEANALYSIS.out.gfastats_purged
    busco_purged        = NFCORE_INTEGRATIONSITEANALYSIS.out.busco_purged
    quast_purged        = NFCORE_INTEGRATIONSITEANALYSIS.out.quast_purged
    purged_assemblies   = NFCORE_INTEGRATIONSITEANALYSIS.out.purged_assemblies
    purgedups_files     = NFCORE_INTEGRATIONSITEANALYSIS.out.purgedups_files
    blast_results       = NFCORE_INTEGRATIONSITEANALYSIS.out.blast_results
    extracted_hits      = NFCORE_INTEGRATIONSITEANALYSIS.out.extracted_hits
    plasmid_gb          = NFCORE_INTEGRATIONSITEANALYSIS.out.plasmid_gb
    plasmid_gff3        = NFCORE_INTEGRATIONSITEANALYSIS.out.plasmid_gff3
    plasmid_features    = NFCORE_INTEGRATIONSITEANALYSIS.out.plasmid_features
    hit_annotations_gff = NFCORE_INTEGRATIONSITEANALYSIS.out.hit_annotations_gff
    hit_annotations_gb  = NFCORE_INTEGRATIONSITEANALYSIS.out.hit_annotations_gb
    mapped_bam          = NFCORE_INTEGRATIONSITEANALYSIS.out.mapped_bam
    mapped_bai          = NFCORE_INTEGRATIONSITEANALYSIS.out.mapped_bai
    coverage_summary    = NFCORE_INTEGRATIONSITEANALYSIS.out.coverage_summary
    igv_report          = NFCORE_INTEGRATIONSITEANALYSIS.out.igv_report
    dotplot_paf         = NFCORE_INTEGRATIONSITEANALYSIS.out.dotplot_paf
    dotplot_gaps        = NFCORE_INTEGRATIONSITEANALYSIS.out.dotplot_gaps
    dotplot_best        = NFCORE_INTEGRATIONSITEANALYSIS.out.dotplot_best
    dotplot_best_fasta  = NFCORE_INTEGRATIONSITEANALYSIS.out.dotplot_best_fasta
    dotplot_genome      = NFCORE_INTEGRATIONSITEANALYSIS.out.dotplot_genome
    dotplot_chromosome  = NFCORE_INTEGRATIONSITEANALYSIS.out.dotplot_chromosome
    dotplot_moddotplot  = NFCORE_INTEGRATIONSITEANALYSIS.out.dotplot_moddotplot
    lifted_annotation   = NFCORE_INTEGRATIONSITEANALYSIS.out.lifted_annotation
    lifted_validation   = NFCORE_INTEGRATIONSITEANALYSIS.out.lifted_validation
    lifted_unmapped     = NFCORE_INTEGRATIONSITEANALYSIS.out.lifted_unmapped
    read_inspector      = NFCORE_INTEGRATIONSITEANALYSIS.out.read_inspector
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PIPELINE OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

output {
    multiqc_report {
        path { _file ->
            "multiqc"
        }
        mode params.publish_dir_mode
    }
    multiqc_data {
        path { _file ->
            "multiqc"
        }
        mode params.publish_dir_mode
    }
    multiqc_plots {
        path { _file ->
            "multiqc"
        }
        mode params.publish_dir_mode
    }
    versions {
        mode params.publish_dir_mode
    }
    read_stats {
        path { meta, _file ->
            "${meta.id}/raw_read_stats"
        }
        mode params.publish_dir_mode
    }
    genomescope {
        path { meta, _file ->
            "${meta.id}/raw_read_stats"
        }
        mode params.publish_dir_mode
    }
    genomescope_plots {
        path { meta, _file ->
            "${meta.id}/raw_read_stats/genomescope2/${meta.id}"
        }
        mode params.publish_dir_mode
    }
    downsample_summary {
        path { meta, _file ->
            "${meta.id}/downsampled_reads_stats"
        }
        mode params.publish_dir_mode
    }
    draft_assemblies {
        path { meta, _primary, _alternate ->
            "${meta.id}/draft_assembly"
        }
        mode params.publish_dir_mode
    }
    asmstats_raw {
        path { meta, _file ->
            "${meta.id}/genome_stats/raw"
        }
        mode params.publish_dir_mode
    }
    gfastats_raw {
        path { meta, _file ->
            "${meta.id}/genome_stats/raw"
        }
        mode params.publish_dir_mode
    }
    busco_raw {
        path { meta, _file ->
            "${meta.id}/genome_stats/raw"
        }
        mode params.publish_dir_mode
    }
    quast_raw {
        path { meta, _file ->
            "${meta.id}/genome_stats/raw"
        }
        mode params.publish_dir_mode
    }
    asmstats_purged {
        path { meta, _file ->
            "${meta.id}/genome_stats/purged"
        }
        mode params.publish_dir_mode
    }
    gfastats_purged {
        path { meta, _file ->
            "${meta.id}/genome_stats/purged"
        }
        mode params.publish_dir_mode
    }
    busco_purged {
        path { meta, _file ->
            "${meta.id}/genome_stats/purged"
        }
        mode params.publish_dir_mode
    }
    quast_purged {
        path { meta, _file ->
            "${meta.id}/genome_stats/purged"
        }
        mode params.publish_dir_mode
    }
    purged_assemblies {
        path { meta, _purged, _pooled_alternate ->
            "${meta.id}/purged_assembly"
        }
        mode params.publish_dir_mode
    }
    purgedups_files {
        path { meta, _file ->
            "${meta.id}/purgedups"
        }
        mode params.publish_dir_mode
    }
    blast_results {
        path { meta, _file ->
            "${meta.id}/integration-results/${meta.db_type}/blast_results"
        }
        mode params.publish_dir_mode
    }
    extracted_hits {
        path { meta, _file ->
            "${meta.id}/integration-results/${meta.db_type}/extracted_hits"
        }
        mode params.publish_dir_mode
    }
    plasmid_gb {
        path { meta, _file ->
            "${meta.id}/plasmid"
        }
        mode params.publish_dir_mode
    }
    plasmid_gff3 {
        path { meta, _file ->
            "${meta.id}/plasmid"
        }
        mode params.publish_dir_mode
    }
    plasmid_features {
        path { meta, _file ->
            "${meta.id}/plasmid"
        }
        mode params.publish_dir_mode
    }
    hit_annotations_gff {
        path { meta, _file ->
            "${meta.id}/integration-results/${meta.db_type}/annotations/${meta.hit_id}"
        }
        mode params.publish_dir_mode
    }
    hit_annotations_gb {
        path { meta, _file ->
            "${meta.id}/integration-results/${meta.db_type}/annotations/${meta.hit_id}"
        }
        mode params.publish_dir_mode
    }
    mapped_bam {
        path { meta, _file ->
            "${meta.id}/integration-results/${meta.db_type}/mapping/${meta.hit_id}"
        }
        mode params.publish_dir_mode
    }
    mapped_bai {
        path { meta, _file ->
            "${meta.id}/integration-results/${meta.db_type}/mapping/${meta.hit_id}"
        }
        mode params.publish_dir_mode
    }
    coverage_summary {
        path { _file ->
            "mapping/alignment_summary"
        }
        mode params.publish_dir_mode
    }
    igv_report {
        path { meta, _file ->
            "${meta.id}/integration-results/${meta.db_type}/mapping/${meta.hit_id}"
        }
        mode params.publish_dir_mode
    }
    dotplot_paf {
        path { meta, _file ->
            "${meta.id}/integration-results/${meta.db_type}/dotplot/${meta.hit_id}"
        }
        mode params.publish_dir_mode
    }
    dotplot_gaps {
        path { meta, _file ->
            "${meta.id}/integration-results/${meta.db_type}/dotplot/${meta.hit_id}"
        }
        mode params.publish_dir_mode
    }
    dotplot_best {
        path { meta, _file ->
            "${meta.id}/integration-results/${meta.db_type}/dotplot/${meta.hit_id}"
        }
        mode params.publish_dir_mode
    }
    dotplot_best_fasta {
        path { meta, _file ->
            "${meta.id}/integration-results/${meta.db_type}/dotplot/${meta.hit_id}"
        }
        mode params.publish_dir_mode
    }
    dotplot_genome {
        path { meta, _file ->
            "${meta.id}/integration-results/${meta.db_type}/dotplot/${meta.hit_id}"
        }
        mode params.publish_dir_mode
    }
    dotplot_chromosome {
        path { meta, _file ->
            "${meta.id}/integration-results/${meta.db_type}/dotplot/${meta.hit_id}"
        }
        mode params.publish_dir_mode
    }
    dotplot_moddotplot {
        path { meta, _file ->
            "${meta.id}/integration-results/${meta.db_type}/dotplot/${meta.hit_id}"
        }
        mode params.publish_dir_mode
    }
    lifted_annotation {
        path { meta, _file ->
            "${meta.id}/integration-results/${meta.db_type}/reference-annotations/${meta.hit_id}"
        }
        mode params.publish_dir_mode
    }
    lifted_validation {
        path { meta, _file ->
            "${meta.id}/integration-results/${meta.db_type}/reference-annotations/${meta.hit_id}"
        }
        mode params.publish_dir_mode
    }
    lifted_unmapped {
        path { meta, _file ->
            "${meta.id}/integration-results/${meta.db_type}/reference-annotations/${meta.hit_id}"
        }
        mode params.publish_dir_mode
    }
    // The emitted value is the read-inspector directory. Publishing into the
    // database folder keeps that directory name, so the session's relative
    // paths still resolve.
    read_inspector {
        path { meta, _dir ->
            "${meta.id}/integration-results/${meta.db_type}"
        }
        mode params.publish_dir_mode
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
