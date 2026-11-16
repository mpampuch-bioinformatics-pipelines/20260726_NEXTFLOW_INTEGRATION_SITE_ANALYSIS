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
    ch_hits   // channel: [ val(meta), path(hit_fastas) ]  (meta carries db_type; hit_fastas may be a list)
    ch_reads  // channel: [ val(id), path(reads) ]         (per-sample reads)
    ch_gff    // channel: [ val(meta), path(gff3) ]        (meta carries hit_id; plasmid annotation)
    ch_loci   // channel: [ val(meta), path(loci_tsv) ]    (alignment vs extract coordinates)
    ch_lifted // channel: [ val(meta), path(gff) ]       (lifted reference genes; may be empty)

    main:
    //
    // One channel element per extracted hit FASTA, tagged with hit_id for unique prefixes.
    // path("all_hits/*.fa") or path("auto_loci_prediction/*.fa") may be a single Path
    // when only one sequence exists — normalize to a list.
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
            ref: [meta, hit]
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
        COVERAGE_STATS.out.coverage.map { _meta, tsv -> tsv }.collect().filter { tsvs -> tsvs instanceof List ? tsvs.size() > 0 : tsvs != null }
    )

    //
    // IGV HTML report per hit that has a plasmid GFF. Mapping does not wait
    // on the GFF or the lifted genes; only this join does. A hit with no
    // plasmid GFF is not reported. Lifted genes are a second track when present.
    //
    ch_gff_keyed = ch_gff.map { meta, gff ->
        ["${meta.id}|${meta.db_type}|${meta.hit_id}", gff]
    }
    ch_lifted_keyed = ch_lifted.map { meta, lifted ->
        ["${meta.id}|${meta.db_type}|${meta.hit_id}", lifted]
    }
    // One table per sample and database, copied onto every hit from that database.
    ch_loci_keyed = ch_loci.map { meta, loci ->
        ["${meta.id}|${meta.db_type}", loci]
    }
    ch_igv = MINIMAP2_MAP.out.bam
        .join(MINIMAP2_MAP.out.index)
        .join(ch_per_hit.ref)
        .map { meta, bam, bai, ref ->
            ["${meta.id}|${meta.db_type}|${meta.hit_id}", meta, bam, bai, ref]
        }
        .join(ch_gff_keyed)
        .map { _key, meta, bam, bai, ref, gff ->
            ["${meta.id}|${meta.db_type}|${meta.hit_id}", meta, bam, bai, ref, gff]
        }
        .join(ch_lifted_keyed, remainder: true)
        .map { _key, meta, bam, bai, ref, gff, lifted ->
            // An empty list is "no lifted track". Nextflow 26 does not accept
            // path(optional: true), and a null from the remainder join is not a file.
            ["${meta.id}|${meta.db_type}", meta, bam, bai, ref, gff, lifted ?: []]
        }
        .combine(ch_loci_keyed, by: 0)
        .map { _key, meta, bam, bai, ref, gff, lifted, loci ->
            // Optional path input: absent lift is an empty list, not null (Nextflow 25.10).
            [meta, bam, bai, ref, gff, lifted ?: [], loci]
        }

    IGV_REPORT(ch_igv)

    // One [meta, png] per snapshot. A single locus emits one path; several
    // BED rows emit a list, and MultiQC needs each file on its own.
    ch_igv_png = IGV_REPORT.out.png.flatMap { meta, pngs ->
        def files = pngs instanceof List ? pngs : [pngs]
        files.collect { png -> [meta, png] }
    }

    emit:
    bam              = MINIMAP2_MAP.out.bam         // channel: [ val(meta), path(bam) ]
    index            = MINIMAP2_MAP.out.index       // channel: [ val(meta), path(bai) ]
    coverage_summary = COVERAGE_SUMMARY.out.summary // channel: path(summary_files)
    igv_report       = IGV_REPORT.out.report        // channel: [ val(meta), path(html) ]
    igv_png          = ch_igv_png                   // channel: [ val(meta), path(png) ]
}
