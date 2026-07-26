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
    versions            = INTEGRATION_SITE_ANALYSIS.out.versions
    read_stats          = INTEGRATION_SITE_ANALYSIS.out.read_stats
    genomescope         = INTEGRATION_SITE_ANALYSIS.out.genomescope
    downsampled_reads   = INTEGRATION_SITE_ANALYSIS.out.downsampled_reads
    downsample_summary  = INTEGRATION_SITE_ANALYSIS.out.downsample_summary
    draft_assemblies    = INTEGRATION_SITE_ANALYSIS.out.draft_assemblies
    asmstats_raw        = INTEGRATION_SITE_ANALYSIS.out.asmstats_raw
    gfastats_raw        = INTEGRATION_SITE_ANALYSIS.out.gfastats_raw
    busco_raw           = INTEGRATION_SITE_ANALYSIS.out.busco_raw
    asmstats_purged     = INTEGRATION_SITE_ANALYSIS.out.asmstats_purged
    gfastats_purged     = INTEGRATION_SITE_ANALYSIS.out.gfastats_purged
    busco_purged        = INTEGRATION_SITE_ANALYSIS.out.busco_purged
    purged_assemblies   = INTEGRATION_SITE_ANALYSIS.out.purged_assemblies
    blast_db            = INTEGRATION_SITE_ANALYSIS.out.blast_db
    blast_results       = INTEGRATION_SITE_ANALYSIS.out.blast_results
    extracted_sequences = INTEGRATION_SITE_ANALYSIS.out.extracted_sequences
    extracted_hits      = INTEGRATION_SITE_ANALYSIS.out.extracted_hits
    mapped_bam          = INTEGRATION_SITE_ANALYSIS.out.mapped_bam
    mapped_bai          = INTEGRATION_SITE_ANALYSIS.out.mapped_bai
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
    versions            = NFCORE_INTEGRATIONSITEANALYSIS.out.versions
    read_stats          = NFCORE_INTEGRATIONSITEANALYSIS.out.read_stats
    genomescope         = NFCORE_INTEGRATIONSITEANALYSIS.out.genomescope
    downsampled_reads   = NFCORE_INTEGRATIONSITEANALYSIS.out.downsampled_reads
    downsample_summary  = NFCORE_INTEGRATIONSITEANALYSIS.out.downsample_summary
    draft_assemblies    = NFCORE_INTEGRATIONSITEANALYSIS.out.draft_assemblies
    asmstats_raw        = NFCORE_INTEGRATIONSITEANALYSIS.out.asmstats_raw
    gfastats_raw        = NFCORE_INTEGRATIONSITEANALYSIS.out.gfastats_raw
    busco_raw           = NFCORE_INTEGRATIONSITEANALYSIS.out.busco_raw
    asmstats_purged     = NFCORE_INTEGRATIONSITEANALYSIS.out.asmstats_purged
    gfastats_purged     = NFCORE_INTEGRATIONSITEANALYSIS.out.gfastats_purged
    busco_purged        = NFCORE_INTEGRATIONSITEANALYSIS.out.busco_purged
    purged_assemblies   = NFCORE_INTEGRATIONSITEANALYSIS.out.purged_assemblies
    blast_db            = NFCORE_INTEGRATIONSITEANALYSIS.out.blast_db
    blast_results       = NFCORE_INTEGRATIONSITEANALYSIS.out.blast_results
    extracted_sequences = NFCORE_INTEGRATIONSITEANALYSIS.out.extracted_sequences
    extracted_hits      = NFCORE_INTEGRATIONSITEANALYSIS.out.extracted_hits
    mapped_bam          = NFCORE_INTEGRATIONSITEANALYSIS.out.mapped_bam
    mapped_bai          = NFCORE_INTEGRATIONSITEANALYSIS.out.mapped_bai
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PIPELINE OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

output {
    multiqc_report {
        mode params.publish_dir_mode
    }
    versions {
        mode params.publish_dir_mode
    }
    read_stats {
        path { meta, _file ->
            "${meta.id}/read_qc"
        }
        mode params.publish_dir_mode
    }
    genomescope {
        path { meta, _file ->
            "${meta.id}/read_qc"
        }
        mode params.publish_dir_mode
    }
    downsampled_reads {
        path { meta, _file ->
            "${meta.id}/downsampled"
        }
        mode params.publish_dir_mode
    }
    downsample_summary {
        path { meta, _file ->
            "${meta.id}/downsampled"
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
    purged_assemblies {
        path { meta, _purged, _haplotigs ->
            "${meta.id}/purged_assembly"
        }
        mode params.publish_dir_mode
    }
    blast_db {
        path { meta, _file ->
            "${meta.id}/${meta.db_type}/blast_db"
        }
        mode params.publish_dir_mode
    }
    blast_results {
        path { meta, _file ->
            "${meta.id}/${meta.db_type}/blast_results"
        }
        mode params.publish_dir_mode
    }
    extracted_sequences {
        path { meta, _file ->
            "${meta.id}/${meta.db_type}/extracted_sequences"
        }
        mode params.publish_dir_mode
    }
    extracted_hits {
        path { meta, _file ->
            "${meta.id}/${meta.db_type}/extracted_hits"
        }
        mode params.publish_dir_mode
    }
    mapped_bam {
        path { meta, _file ->
            "${meta.id}/${meta.db_type}/mapping"
        }
        mode params.publish_dir_mode
    }
    mapped_bai {
        path { meta, _file ->
            "${meta.id}/${meta.db_type}/mapping"
        }
        mode params.publish_dir_mode
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
