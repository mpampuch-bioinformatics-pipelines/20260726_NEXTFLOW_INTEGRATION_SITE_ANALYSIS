//
// Dot plots of each extracted hit against its reference.
// Whole-reference mummerplot, plus a second plot of the reference contig
// with the most minimap2 alignment-block bases. modDotPlot is the identity check.
//

include { MINIMAP2_ALIGN as MINIMAP2_DOTPLOT } from '../../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_REFERENCE } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_QUERY     } from '../../../modules/nf-core/samtools/faidx/main'
include { PAF_HIT_REFERENCE                          } from '../../../modules/local/paf_hit_reference/main'
include { MUMMER_DOTPLOT                             } from '../../../modules/local/mummer_dotplot/main'
include { MODDOTPLOT                                 } from '../../../modules/local/moddotplot/main'

workflow HIT_REFERENCE_DOTPLOT {
    take:
    ch_hits // channel: [ val(meta), path(hit_fastas) ]  (hit_fastas may be a list; meta carries db_type)
    ch_reference // channel: [ val(meta), path(reference) ]

    main:
    // One element per extracted FASTA. A single path is not a list.
    ch_per_hit = ch_hits
        .map { meta, hits ->
            def hit_list = hits instanceof List ? hits : [hits]
            [meta, hit_list]
        }
        .transpose()
        .map { meta, hit ->
            [meta.id, meta + [hit_id: hit.baseName], hit]
        }

    ch_ref = ch_reference.map { meta, reference ->
        [meta.id, reference]
    }

    ch_paired = ch_per_hit
        .combine(ch_ref, by: 0)
        .multiMap { _id, meta, hit, reference ->
            query: [meta, hit]
            reference: [meta, reference]
        }

    // PAF with CIGAR. Preset is ext.args on MINIMAP2_DOTPLOT (-x asm5 --secondary=no).
    MINIMAP2_DOTPLOT(ch_paired.query, ch_paired.reference, false, 'bai', true, false)

    // Indexes sit beside the FASTAs that modDotPlot reads. No prior index.
    SAMTOOLS_FAIDX_REFERENCE(
        ch_paired.reference.map { meta, reference -> [meta, reference, []] },
        false
    )
    SAMTOOLS_FAIDX_QUERY(
        ch_paired.query.map { meta, query -> [meta, query, []] },
        false
    )

    PAF_HIT_REFERENCE(
        MINIMAP2_DOTPLOT.out.paf.join(ch_paired.reference)
    )

    MUMMER_DOTPLOT(
        ch_paired.reference
            .join(ch_paired.query)
            .join(PAF_HIT_REFERENCE.out.best_fasta)
    )

    MODDOTPLOT(
        PAF_HIT_REFERENCE.out.pairs
            .join(ch_paired.reference)
            .join(SAMTOOLS_FAIDX_REFERENCE.out.fai)
            .join(ch_paired.query)
            .join(SAMTOOLS_FAIDX_QUERY.out.fai)
    )

    emit:
    paf              = MINIMAP2_DOTPLOT.out.paf
    gaps             = PAF_HIT_REFERENCE.out.gaps
    best_chromosome  = PAF_HIT_REFERENCE.out.best_chromosome
    best_fasta       = PAF_HIT_REFERENCE.out.best_fasta
    genome_png       = MUMMER_DOTPLOT.out.genome_png
    chromosome_png   = MUMMER_DOTPLOT.out.chromosome_png
    moddotplot       = MODDOTPLOT.out.moddotplot
}
