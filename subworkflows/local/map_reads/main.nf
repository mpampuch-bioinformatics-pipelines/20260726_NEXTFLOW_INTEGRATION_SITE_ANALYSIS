//
// Map reads to each extracted integration-site hit (nxf-minimapper flow):
// minimap2 BAM, Medaka consensus/variants, coverage summary, IGV reports.
//

include { MINIMAP2_ALIGN as MINIMAP2_MAP } from '../../../modules/nf-core/minimap2/align/main'
include { COVERAGE_STATS                 } from '../../../modules/local/coverage_stats/main'
include { COVERAGE_SUMMARY               } from '../../../modules/local/coverage_summary/main'
include { IGV_REPORT                     } from '../../../modules/local/igv_report/main'

workflow MAP_READS {

    take:
    ch_hits  // channel: [ val(meta), path(hit_fastas) ]  (meta carries db_type; hit_fastas may be a list)
    ch_reads // channel: [ val(id), path(reads) ]         (per-sample reads)

    main:
    //
    // One channel element per extracted hit FASTA, tagged with hit_id for unique prefixes.
    // path("hits/*.fa") may be a single Path when only one hit exists — normalize to a list.
    //
    ch_per_hit = ch_hits
        .map { meta, hits ->
            def hit_list = hits instanceof List ? hits : [hits]
            [meta, hit_list]
        }
        .transpose()
        .map { meta, hit ->
            def hit_id = hit.baseName
            [meta.id, meta + [hit_id: hit_id], hit]
        }
        .combine(ch_reads, by: 0)
        .multiMap { _id, meta, hit, reads ->
            reads: [meta, reads]
            ref:   [meta, hit]
            pair:  [meta, hit, reads]
        }

    //
    // Parallel: map reads to hit, call variants, build consensus (scaffold pattern).
    //
    MINIMAP2_MAP(ch_per_hit.reads, ch_per_hit.ref, true, 'bai', false, false)

    //
    // Coverage + consensus quality sparklines, then one HTML summary across all hits.
    //
    COVERAGE_STATS(
        MINIMAP2_MAP.out.bam.join(MINIMAP2_MAP.out.index)
    )

    // Skip summary when there were no hits (empty collect would leave R with no inputs).
    COVERAGE_SUMMARY(
        COVERAGE_STATS.out.coverage
            .map { _meta, tsv -> tsv }
            .collect()
            .filter { tsvs -> tsvs instanceof List ? tsvs.size() > 0 : tsvs != null }
    )

    //
    // IGV HTML report per hit (optional VCF track controlled by params.variants).
    //
    ch_igv = MINIMAP2_MAP.out.bam
        .join(MINIMAP2_MAP.out.index)
        .join(ch_per_hit.ref)

    IGV_REPORT(ch_igv)

    emit:
    bam              = MINIMAP2_MAP.out.bam           // channel: [ val(meta), path(bam) ]
    index            = MINIMAP2_MAP.out.index         // channel: [ val(meta), path(bai) ]
    coverage_summary = COVERAGE_SUMMARY.out.summary   // channel: path(summary_files)
    igv_report       = IGV_REPORT.out.report          // channel: [ val(meta), path(html) ]
}
