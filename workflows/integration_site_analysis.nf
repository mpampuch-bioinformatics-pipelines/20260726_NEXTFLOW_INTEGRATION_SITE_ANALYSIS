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
include { ASSEMBLY_QC as ASSEMBLY_QC_RAW                } from '../subworkflows/local/assembly_qc'
include { ASSEMBLY_QC as ASSEMBLY_QC_PURGED             } from '../subworkflows/local/assembly_qc'
include { PURGE_ASSEMBLY                                } from '../subworkflows/local/purge_assembly'
include { BLAST_DATABASES                               } from '../subworkflows/local/blast_databases'
include { BLAST_SEARCH                                  } from '../subworkflows/local/blast_search'
include { EXTRACT_SEQUENCES                             } from '../subworkflows/local/extract_sequences'
include { MAP_READS                                     } from '../subworkflows/local/map_reads'
include { READ_INSPECTOR                                } from '../subworkflows/local/read_inspector/main'
include { JUNCTION_QC                                   } from '../subworkflows/local/junction_qc/main'
include { HIT_REFERENCE_DOTPLOT                       } from '../subworkflows/local/hit_reference_dotplot/main'
include { LIFT_REFERENCE_ANNOTATIONS                  } from '../subworkflows/local/lift_reference_annotations/main'
include { PLASMID_ANNOTATION                            } from '../subworkflows/local/plasmid_annotation'
include { CONCAT_FASTA as CONCAT_QUERIES                } from '../modules/local/concat_fasta/main'
include { GUNZIP                                        } from '../modules/nf-core/gunzip/main'
include { SEQKIT_STATS as REFERENCE_SIZE                } from '../modules/nf-core/seqkit/stats/main'
include { BUSCO_DOWNLOAD                                } from '../modules/nf-core/busco/download/main'

// Shared processes are keyed by absolute path, not sample id. The task hash
// includes meta.id, so two samples that point at one file would otherwise
// run the process twice.
def sharedMeta(key, label) {
    def hash = java.security.MessageDigest
        .getInstance('MD5')
        .digest(key.getBytes('UTF-8'))
        .encodeHex()
        .toString()
        .substring(0, 8)

    def base = key.contains('|') ? label : new File(key).name
    def safe = base.replaceAll(/[^A-Za-z0-9._-]/, '_')

    [id: "${safe}_${hash}", key: key]
}

def fanShared(shared_files, key_to_sample) {
    shared_files
        .map { meta, file -> [meta.key, file] }
        .combine(key_to_sample, by: 0)
        .map { _key, file, sample_meta -> [sample_meta, file] }
}

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
    // Named BUSCO lineages are downloaded once. Auto-lineage stays inside BUSCO.
    //
    ch_named_lineages = ch_samplesheet
        .map { meta, _reads -> meta.busco_lineage ?: params.busco_lineage ?: 'auto_euk' }
        .unique()
        .filter { lineage -> !(lineage in ['auto', 'auto_prok', 'auto_euk']) }
    BUSCO_DOWNLOAD(ch_named_lineages)
    // [lineage, busco_downloads]. The name is part of the process output, so a
    // later download cannot be paired with an earlier request.
    ch_lineage_db = BUSCO_DOWNLOAD.out.download_dir

    //
    // SUBWORKFLOW: QC the input reads once per distinct FASTQ.
    //
    ch_reads_keyed = ch_samplesheet.map { meta, reads ->
        def reads_file = file(reads, checkIfExists: true)
        [reads_file.toAbsolutePath().toString(), meta, reads_file]
    }
    READ_QC(
        ch_reads_keyed.map { key, _meta, reads_file -> [sharedMeta(key, 'reads'), reads_file] }.unique { shared, _reads_file -> shared.key }
    )
    ch_multiqc_files = ch_multiqc_files.mix(READ_QC.out.multiqc_files)
    ch_reads_by_key = ch_reads_keyed.map { key, meta, _reads_file -> [key, meta] }
    ch_read_stats = fanShared(READ_QC.out.read_stats, ch_reads_by_key)
    ch_genomescope = fanShared(READ_QC.out.genomescope, ch_reads_by_key)
    ch_genomescope_plots = fanShared(READ_QC.out.genomescope_plots, ch_reads_by_key)

    //
    // Reference length is computed once per reference file (nf-core/seqkit/stats).
    //
    ch_reference_keyed = ch_samplesheet.map { meta, _reads ->
        def reference = file(meta.reference ?: params.reference, checkIfExists: true)
        [reference.toAbsolutePath().toString(), meta, reference]
    }
    REFERENCE_SIZE(
        ch_reference_keyed.map { key, _meta, reference -> [sharedMeta(key, 'reference'), reference] }.unique { shared, _reference -> shared.key }
    )
    ch_genome_size = REFERENCE_SIZE.out.stats.map { shared, tsv -> [shared.key, tsv] }
    ch_reference = ch_reference_keyed.map { _key, meta, reference -> [meta, reference] }

    //
    // SUBWORKFLOW: Downsample once per reads file, reference, coverage, and seed.
    //
    ch_downsample_keyed = ch_samplesheet
        .map { meta, reads ->
            def reference = file(meta.reference ?: params.reference, checkIfExists: true)
            def reads_file = file(reads, checkIfExists: true)
            def ref_key = reference.toAbsolutePath().toString()
            def coverage = meta.target_coverage ?: params.target_coverage
            def job_key = "${reads_file.toAbsolutePath()}|${ref_key}|${coverage}|${params.seqkit_seed}"
            [ref_key, job_key, meta, reads_file]
        }
        .combine(ch_genome_size, by: 0)
        .map { _ref_key, job_key, meta, reads_file, size_tsv ->
            [job_key, meta, reads_file, size_tsv]
        }
    DOWNSAMPLE_READS(
        ch_downsample_keyed.map { job_key, meta, reads_file, size_tsv ->
            def coverage = meta.target_coverage ?: params.target_coverage
            [sharedMeta(job_key, 'downsample') + [target_coverage: coverage], reads_file, size_tsv]
        }.unique { shared, _reads_file, _size_tsv -> shared.key }
    )
    ch_downsample_by_key = ch_downsample_keyed.map { job_key, meta, _reads_file, _size_tsv -> [job_key, meta] }
    ch_downsampled = fanShared(DOWNSAMPLE_READS.out.reads, ch_downsample_by_key)
    ch_downsample_summary = fanShared(DOWNSAMPLE_READS.out.summary, ch_downsample_by_key)

    //
    // Samples with both primary and alternate skip hifiasm. Gzipped user FASTAs
    // go through nf-core/gunzip; uncompressed files pass through unchanged.
    //
    ch_assembly_input = ch_downsampled.branch { meta, _reads ->
        assemble: !meta.primary
        provided: meta.primary
    }

    // One hifiasm run per distinct downsampled FASTQ. Purge still branches per sample.
    ch_assemble_keyed = ch_assembly_input.assemble.map { meta, reads ->
        def key = file(reads).toAbsolutePath().toString()
        [key, meta, reads]
    }
    DRAFT_ASSEMBLY(
        ch_assemble_keyed.map { key, _meta, reads -> [sharedMeta(key, 'assembly'), reads] }.unique { shared, _reads -> shared.key }
    )
    ch_draft_from_hifiasm = DRAFT_ASSEMBLY.out.assemblies
        .map { shared, primary, alternate -> [shared.key, primary, alternate] }
        .combine(ch_assemble_keyed.map { key, meta, _reads -> [key, meta] }, by: 0)
        .map { _key, primary, alternate, meta -> [meta, primary, alternate] }

    ch_user_fastas = ch_assembly_input.provided
        .flatMap { meta, _reads ->
            [[meta + [assembly: 'primary'], file(meta.primary, checkIfExists: true)], [meta + [assembly: 'alternate'], file(meta.alternate, checkIfExists: true)]]
        }
        .branch { _meta, fasta ->
            gzipped: fasta.name.endsWith('.gz')
            plain: true
        }

    GUNZIP(ch_user_fastas.gzipped)

    ch_user_unzipped = GUNZIP.out.gunzip.mix(ch_user_fastas.plain)

    ch_user_primary = ch_user_unzipped
        .filter { meta, _fasta -> meta.assembly == 'primary' }
        .map { meta, fasta ->
            def clean = meta.clone()
            clean.remove('assembly')
            [meta.id, clean, fasta]
        }
    ch_user_alternate = ch_user_unzipped
        .filter { meta, _fasta -> meta.assembly == 'alternate' }
        .map { meta, fasta -> [meta.id, fasta] }

    ch_user_assemblies = ch_user_primary
        .join(ch_user_alternate)
        .map { _id, meta, primary, alternate -> [meta, primary, alternate] }

    ch_draft_assemblies = ch_draft_from_hifiasm.mix(ch_user_assemblies)

    //
    // SUBWORKFLOW: QC the raw draft assemblies
    //
    GENOME_STATISTICS_RAW(ch_draft_assemblies, 'raw', ch_lineage_db)
    ch_multiqc_files = ch_multiqc_files.mix(GENOME_STATISTICS_RAW.out.multiqc_files)

    ASSEMBLY_QC_RAW(ch_draft_assemblies, ch_reference, 'raw')
    ch_multiqc_files = ch_multiqc_files.mix(ASSEMBLY_QC_RAW.out.multiqc_files)

    //
    // Decide per-sample whether to purge (samplesheet override or global param).
    //
    ch_draft_branched = ch_draft_assemblies.branch { meta, _primary, _alternate ->
        purge: (meta.purge != null ? meta.purge : params.run_purging)
        keep: !(meta.purge != null ? meta.purge : params.run_purging)
    }

    //
    // SUBWORKFLOW: Optional purge_dups on the primary assembly, then re-QC
    //
    PURGE_ASSEMBLY(ch_draft_branched.purge, ch_downsampled)
    GENOME_STATISTICS_PURGED(PURGE_ASSEMBLY.out.assemblies, 'purged', ch_lineage_db)
    ch_multiqc_files = ch_multiqc_files.mix(GENOME_STATISTICS_PURGED.out.multiqc_files)

    ASSEMBLY_QC_PURGED(PURGE_ASSEMBLY.out.assemblies, ch_reference, 'purged')
    ch_multiqc_files = ch_multiqc_files.mix(ASSEMBLY_QC_PURGED.out.multiqc_files)

    // Assemblies feeding BLAST: purged primary + (original alternate ∪ haplotigs), or kept pair (p_ctg + a_ctg).
    ch_for_blast = PURGE_ASSEMBLY.out.assemblies.mix(ch_draft_branched.keep)

    //
    // SUBWORKFLOW: Build BLAST databases (primary / alternate / combined)
    //
    BLAST_DATABASES(ch_for_blast)

    //
    // Resolve and merge the construct query FASTAs per sample.
    //
    ch_query_keyed = ch_samplesheet.map { meta, _reads ->
        def queries = meta.blast_queries ?: params.blast_queries
        def query_path = file(queries, checkIfExists: true)
        def key = query_path.toAbsolutePath().toString()
        def query_files = query_path.isDirectory() ? files("${queries}/*.{fa,fasta,fna}") : [query_path]
        [key, meta.id, query_files]
    }
    CONCAT_QUERIES(
        ch_query_keyed.map { key, _id, query_files -> [sharedMeta(key, 'queries'), query_files] }.unique { shared, _query_files -> shared.key }
    )
    ch_query = CONCAT_QUERIES.out.fasta
        .map { shared, fasta -> [shared.key, fasta] }
        .combine(ch_query_keyed.map { key, id, _query_files -> [key, id] }, by: 0)
        .map { _key, fasta, id -> [id, fasta] }

    //
    // SUBWORKFLOW: BLAST the queries against each database
    //
    BLAST_SEARCH(BLAST_DATABASES.out.db, ch_query)

    //
    // SUBWORKFLOW: Extract hit sequences with flanking context
    //
    EXTRACT_SEQUENCES(BLAST_SEARCH.out.results, BLAST_DATABASES.out.db)

    // auto_loci is one window per insertion. all_hits keeps every query feature.
    ch_downstream_hits = params.downstream_hits == 'all_hits'
        ? EXTRACT_SEQUENCES.out.all_hits
        : EXTRACT_SEQUENCES.out.auto_loci

    //
    // SUBWORKFLOW: Optional plasmid feature annotation of the selected hit set
    //
    ch_plasmid = ch_samplesheet
        .map { meta, _reads ->
            def plasmid = meta.plasmid ?: params.plasmid
            plasmid ? [meta, file(plasmid, checkIfExists: true)] : null
        }
        .filter { it != null }

    PLASMID_ANNOTATION(ch_plasmid, ch_downstream_hits)

    //
    // SUBWORKFLOW: Lift reference genes onto each selected hit when an
    // annotation file is set. Samples without one are not lifted.
    //
    ch_annotation = ch_samplesheet
        .map { meta, _reads ->
            def annotation = meta.annotation ?: params.annotation
            annotation ? [meta.id, file(annotation, checkIfExists: true)] : null
        }
        .filter { it != null }

    LIFT_REFERENCE_ANNOTATIONS(ch_downstream_hits, ch_reference, ch_annotation)

    //
    // SUBWORKFLOW: Map reads to each selected hit (minimap2 + Medaka + coverage + IGV)
    //
    ch_reads_by_id = ch_samplesheet.map { meta, reads -> [meta.id, reads] }
    // auto_loci rows are the clustered insertion. all_hits rows are one BLAST alignment each.
    ch_locus_table = params.downstream_hits == 'all_hits'
        ? EXTRACT_SEQUENCES.out.hit_intervals
        : EXTRACT_SEQUENCES.out.loci_summary
    MAP_READS(
        ch_downstream_hits,
        ch_reads_by_id,
        PLASMID_ANNOTATION.out.hit_gff3,
        ch_locus_table,
        LIFT_REFERENCE_ANNOTATIONS.out.lifted,
    )

    //
    // SUBWORKFLOW: One IGV/JBrowse directory per sample and assembly.
    // Uses the original FASTQ. The per-hit HTML report above is separate.
    //
    READ_INSPECTOR(ch_downstream_hits, ch_reads_by_id)

    //
    // Junction metrics from that BAM and the locus table. The table does not
    // decide trust versus revisit.
    //
    JUNCTION_QC(READ_INSPECTOR.out.inspector, ch_locus_table)

    //
    // SUBWORKFLOW: Dot plots of each selected hit against the reference
    //
    HIT_REFERENCE_DOTPLOT(ch_downstream_hits, ch_reference)

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
    read_stats          = ch_read_stats // channel: [ val(meta), path(tsv) ]
    genomescope         = ch_genomescope // channel: [ val(meta), path(summary) ]
    genomescope_plots   = ch_genomescope_plots // channel: [ val(meta), path(file) ]
    downsample_summary  = ch_downsample_summary // channel: [ val(meta), path(tsv) ]
    draft_assemblies    = ch_draft_assemblies // channel: [ val(meta), path(primary), path(alternate) ]
    asmstats_raw        = GENOME_STATISTICS_RAW.out.asmstats // channel: [ val(meta), path(stats) ]
    gfastats_raw        = GENOME_STATISTICS_RAW.out.gfastats // channel: [ val(meta), path(summary) ]
    busco_raw           = GENOME_STATISTICS_RAW.out.busco // channel: [ val(meta), path(summary) ]
    quast_raw           = ASSEMBLY_QC_RAW.out.results // channel: [ val(meta), path(dir) ]
    asmstats_purged     = GENOME_STATISTICS_PURGED.out.asmstats // channel: [ val(meta), path(stats) ]
    gfastats_purged     = GENOME_STATISTICS_PURGED.out.gfastats // channel: [ val(meta), path(summary) ]
    busco_purged        = GENOME_STATISTICS_PURGED.out.busco // channel: [ val(meta), path(summary) ]
    quast_purged        = ASSEMBLY_QC_PURGED.out.results // channel: [ val(meta), path(dir) ]
    purged_assemblies   = PURGE_ASSEMBLY.out.assemblies // channel: [ val(meta), path(purged), path(pooled_alternate) ]
    purgedups_files     = PURGE_ASSEMBLY.out.purgedups_files // channel: [ val(meta), path(file) ]
    blast_results       = BLAST_SEARCH.out.results // channel: [ val(meta), path(results) ]
    extracted_hits      = EXTRACT_SEQUENCES.out.all_hits.mix(EXTRACT_SEQUENCES.out.auto_loci).mix(EXTRACT_SEQUENCES.out.loci_summary) // channel: [ val(meta), path(hit_or_locus_files) ]
    plasmid_gb          = PLASMID_ANNOTATION.out.plasmid_gb // channel: [ val(meta), path(gb) ]
    plasmid_gff3        = PLASMID_ANNOTATION.out.plasmid_gff3 // channel: [ val(meta), path(gff3) ]
    plasmid_features    = PLASMID_ANNOTATION.out.plasmid_features_fasta // channel: [ val(meta), path(fasta) ]
    hit_annotations_gff = PLASMID_ANNOTATION.out.hit_gff3 // channel: [ val(meta), path(gff3) ]
    hit_annotations_gb  = PLASMID_ANNOTATION.out.hit_gb // channel: [ val(meta), path(gb) ]
    mapped_bam          = MAP_READS.out.bam // channel: [ val(meta), path(bam) ]
    mapped_bai          = MAP_READS.out.index // channel: [ val(meta), path(bai) ]
    coverage_summary    = MAP_READS.out.coverage_summary // channel: path(summary_files)
    igv_report          = MAP_READS.out.igv_report // channel: [ val(meta), path(html) ]
    dotplot_paf         = HIT_REFERENCE_DOTPLOT.out.paf // channel: [ val(meta), path(paf) ]
    dotplot_gaps        = HIT_REFERENCE_DOTPLOT.out.gaps // channel: [ val(meta), path(tsv) ]
    dotplot_best        = HIT_REFERENCE_DOTPLOT.out.best_chromosome // channel: [ val(meta), path(txt) ]
    dotplot_best_fasta  = HIT_REFERENCE_DOTPLOT.out.best_fasta // channel: [ val(meta), path(fasta) ]
    dotplot_genome      = HIT_REFERENCE_DOTPLOT.out.genome_png // channel: [ val(meta), path(png) ]
    dotplot_chromosome  = HIT_REFERENCE_DOTPLOT.out.chromosome_png // channel: [ val(meta), path(png) ]
    dotplot_moddotplot  = HIT_REFERENCE_DOTPLOT.out.moddotplot // channel: [ val(meta), path(dir) ]
    lifted_annotation   = LIFT_REFERENCE_ANNOTATIONS.out.lifted // channel: [ val(meta), path(gff) ]
    lifted_validation   = LIFT_REFERENCE_ANNOTATIONS.out.validation // channel: [ val(meta), path(tsv) ]
    lifted_unmapped     = LIFT_REFERENCE_ANNOTATIONS.out.unmapped // channel: [ val(meta), path(txt) ]
    read_inspector      = READ_INSPECTOR.out.inspector // channel: [ val(meta), path(dir) ]
    junction_qc         = JUNCTION_QC.out.metrics // channel: [ val(meta), path(tsv) ]
}
