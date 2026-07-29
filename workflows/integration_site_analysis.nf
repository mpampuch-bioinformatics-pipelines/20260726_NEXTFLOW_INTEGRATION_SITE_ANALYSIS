/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { READ_QC                                       } from '../subworkflows/local/read_qc'
include { DOWNSAMPLE_READS                              } from '../subworkflows/local/downsample_reads'
include { DRAFT_ASSEMBLY                                } from '../subworkflows/local/draft_assembly'
include { GENOME_STATISTICS as GENOME_STATISTICS_RAW    } from '../subworkflows/local/genome_statistics'
include { GENOME_STATISTICS as GENOME_STATISTICS_PURGED } from '../subworkflows/local/genome_statistics'
include { PURGE_ASSEMBLY                                } from '../subworkflows/local/purge_assembly'
include { BLAST_DATABASES                               } from '../subworkflows/local/blast_databases'
include { BLAST_SEARCH                                  } from '../subworkflows/local/blast_search'
include { EXTRACT_SEQUENCES                             } from '../subworkflows/local/extract_sequences'
include { MAP_READS                                     } from '../subworkflows/local/map_reads'
include { PLASMID_ANNOTATION                            } from '../subworkflows/local/plasmid_annotation'
include { CONCAT_FASTA as CONCAT_QUERIES                } from '../modules/local/concat_fasta/main'

include { MULTIQC                                       } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                              } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                        } from '../subworkflows/local/utils_nfcore_pipeline_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow INTEGRATION_SITE_ANALYSIS {
    take:
    ch_samplesheet // channel: [ val(meta), path(reads) ]

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // Resolve the per-sample reference genome (samplesheet override or global param).
    //
    ch_reads_reference = ch_samplesheet.map { meta, reads ->
        def reference = meta.reference ?: params.reference
        [meta, reads, file(reference, checkIfExists: true)]
    }

    //
    // SUBWORKFLOW: QC the input reads (kmer profile + read stats)
    //
    READ_QC(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(READ_QC.out.multiqc_files)

    //
    // SUBWORKFLOW: Downsample reads to ~target_coverage based on reference size
    //
    DOWNSAMPLE_READS(ch_reads_reference)
    ch_downsampled = DOWNSAMPLE_READS.out.reads

    //
    // SUBWORKFLOW: Draft assembly with hifiasm
    //
    DRAFT_ASSEMBLY(ch_downsampled)

    //
    // SUBWORKFLOW: QC the raw draft assemblies
    //
    GENOME_STATISTICS_RAW(DRAFT_ASSEMBLY.out.assemblies, 'raw')
    ch_multiqc_files = ch_multiqc_files.mix(GENOME_STATISTICS_RAW.out.multiqc_files)

    //
    // Decide per-sample whether to purge (samplesheet override or global param).
    //
    ch_draft_branched = DRAFT_ASSEMBLY.out.assemblies.branch { meta, _primary, _alternate ->
        purge: (meta.purge != null ? meta.purge : params.run_purging)
        keep: !(meta.purge != null ? meta.purge : params.run_purging)
    }

    //
    // SUBWORKFLOW: Optional purge_dups on the primary assembly, then re-QC
    //
    PURGE_ASSEMBLY(ch_draft_branched.purge, ch_downsampled)
    GENOME_STATISTICS_PURGED(PURGE_ASSEMBLY.out.assemblies, 'purged')
    ch_multiqc_files = ch_multiqc_files.mix(GENOME_STATISTICS_PURGED.out.multiqc_files)

    // Assemblies feeding BLAST: purged pair (purged + haplotigs) or kept pair (p_ctg + a_ctg).
    ch_for_blast = PURGE_ASSEMBLY.out.assemblies.mix(ch_draft_branched.keep)

    //
    // SUBWORKFLOW: Build BLAST databases (primary / alternate / combined)
    //
    BLAST_DATABASES(ch_for_blast)

    //
    // Resolve and merge the construct query FASTAs per sample.
    //
    ch_query_input = ch_samplesheet.map { meta, _reads ->
        def queries = meta.blast_queries ?: params.blast_queries
        def query_path = file(queries)
        def query_files = query_path.isDirectory() ? files("${queries}/*.{fa,fasta,fna}") : [query_path]
        [[id: meta.id], query_files]
    }
    CONCAT_QUERIES(ch_query_input)
    ch_query = CONCAT_QUERIES.out.fasta.map { meta, fasta -> [meta.id, fasta] }

    //
    // SUBWORKFLOW: BLAST the queries against each database
    //
    BLAST_SEARCH(BLAST_DATABASES.out.db, ch_query)

    //
    // SUBWORKFLOW: Extract hit sequences with flanking context
    //
    EXTRACT_SEQUENCES(BLAST_SEARCH.out.results, BLAST_DATABASES.out.db)

    //
    // SUBWORKFLOW: Optional plasmid feature annotation of extracted hits
    //
    ch_plasmid = ch_samplesheet.map { meta, _reads ->
        def plasmid = meta.plasmid ?: params.plasmid
        plasmid ? [[id: meta.id], file(plasmid, checkIfExists: true)] : null
    }.filter { it != null }

    PLASMID_ANNOTATION(ch_plasmid, EXTRACT_SEQUENCES.out.hits)

    //
    // SUBWORKFLOW: Map reads to each extracted hit (minimap2 + Medaka + coverage + IGV)
    //
    ch_reads_by_id = ch_samplesheet.map { meta, reads -> [meta.id, reads] }
    MAP_READS(EXTRACT_SEQUENCES.out.hits, ch_reads_by_id)

    //
    // Collate and save software versions
    //
    def topic_versions = Channel
        .topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'pipeline_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? channel.fromPath(params.multiqc_config, checkIfExists: true)
        : channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : channel.empty()

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report      = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    multiqc_data        = MULTIQC.out.data // channel: path(multiqc_data)
    multiqc_plots       = MULTIQC.out.plots // channel: path(multiqc_plots)
    versions            = ch_versions // channel: [ path(versions.yml) ]
    read_stats          = READ_QC.out.read_stats // channel: [ val(meta), path(tsv) ]
    genomescope         = READ_QC.out.genomescope // channel: [ val(meta), path(summary) ]
    genomescope_plots   = READ_QC.out.genomescope_plots // channel: [ val(meta), path(file) ]
    downsample_summary  = DOWNSAMPLE_READS.out.summary // channel: [ val(meta), path(tsv) ]
    draft_assemblies    = DRAFT_ASSEMBLY.out.assemblies // channel: [ val(meta), path(primary), path(alternate) ]
    asmstats_raw        = GENOME_STATISTICS_RAW.out.asmstats // channel: [ val(meta), path(stats) ]
    gfastats_raw        = GENOME_STATISTICS_RAW.out.gfastats // channel: [ val(meta), path(summary) ]
    busco_raw           = GENOME_STATISTICS_RAW.out.busco // channel: [ val(meta), path(summary) ]
    asmstats_purged     = GENOME_STATISTICS_PURGED.out.asmstats // channel: [ val(meta), path(stats) ]
    gfastats_purged     = GENOME_STATISTICS_PURGED.out.gfastats // channel: [ val(meta), path(summary) ]
    busco_purged        = GENOME_STATISTICS_PURGED.out.busco // channel: [ val(meta), path(summary) ]
    purged_assemblies   = PURGE_ASSEMBLY.out.assemblies // channel: [ val(meta), path(purged), path(haplotigs) ]
    purgedups_files     = PURGE_ASSEMBLY.out.purgedups_files // channel: [ val(meta), path(file) ]
    blast_results       = BLAST_SEARCH.out.results // channel: [ val(meta), path(results) ]
    extracted_hits      = EXTRACT_SEQUENCES.out.hits // channel: [ val(meta), path(per_hit_fastas) ]
    plasmid_gb          = PLASMID_ANNOTATION.out.plasmid_gb // channel: [ val(meta), path(gb) ]
    plasmid_gff3        = PLASMID_ANNOTATION.out.plasmid_gff3 // channel: [ val(meta), path(gff3) ]
    plasmid_features    = PLASMID_ANNOTATION.out.plasmid_features_fasta // channel: [ val(meta), path(fasta) ]
    hit_annotations_gff = PLASMID_ANNOTATION.out.hit_gff3 // channel: [ val(meta), path(gff3) ]
    hit_annotations_gb  = PLASMID_ANNOTATION.out.hit_gb // channel: [ val(meta), path(gb) ]
    mapped_bam          = MAP_READS.out.bam // channel: [ val(meta), path(bam) ]
    mapped_bai          = MAP_READS.out.index // channel: [ val(meta), path(bai) ]
    coverage_summary    = MAP_READS.out.coverage_summary // channel: path(summary_files)
    igv_report          = MAP_READS.out.igv_report // channel: [ val(meta), path(html) ]
}
