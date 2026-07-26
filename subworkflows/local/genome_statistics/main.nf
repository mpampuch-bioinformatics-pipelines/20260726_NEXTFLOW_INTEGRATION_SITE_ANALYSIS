//
// Assembly QC statistics: asmstats, gfastats and BUSCO
//

include { ASMSTATS    } from '../../../modules/sanger-tol/asmstats/main'
include { GFASTATS    } from '../../../modules/nf-core/gfastats/main'
include { BUSCO_BUSCO } from '../../../modules/nf-core/busco/busco/main'

workflow GENOME_STATISTICS {

    take:
    ch_assemblies // channel: [ val(meta), path(primary), path(alternate) ]
    val_stage     // string: assembly stage label (e.g. 'raw' or 'purged')

    main:
    ch_multiqc_files = channel.empty()

    //
    // Evaluate the primary and alternate assemblies as independent haplotypes.
    // `stage` and `assembly` are stamped into the meta so downstream file names
    // and publish paths stay unique across the raw/purged QC passes.
    //
    ch_split = ch_assemblies
        .flatMap { meta, primary, alternate ->
            def meta_stage = meta + [stage: val_stage]
            def out = [[meta_stage + [assembly: 'primary'], primary]]
            if (alternate) {
                out << [meta_stage + [assembly: 'alternate'], alternate]
            }
            return out
        }
        .filter { _meta, assembly -> assembly }

    //
    // MODULE: contiguity statistics
    //
    ASMSTATS(ch_split)
    ch_multiqc_files = ch_multiqc_files.mix(ASMSTATS.out.stats.map { _meta, stats -> stats })

    GFASTATS(
        ch_split,
        'fasta',
        '',
        '',
        [[:], []],
        [[:], []],
        [[:], []],
        [[:], []],
    )
    ch_multiqc_files = ch_multiqc_files.mix(GFASTATS.out.assembly_summary.map { _meta, summary -> summary })

    //
    // MODULE: gene-completeness assessment with BUSCO
    //
    def busco_lineage = params.busco_lineage ?: 'auto_euk'
    BUSCO_BUSCO(
        ch_split,
        'genome',
        busco_lineage,
        [],
        [],
        true,
    )
    ch_multiqc_files = ch_multiqc_files.mix(BUSCO_BUSCO.out.short_summaries_txt.map { _meta, summary -> summary })

    emit:
    asmstats      = ASMSTATS.out.stats               // channel: [ val(meta), path(stats) ]
    gfastats      = GFASTATS.out.assembly_summary    // channel: [ val(meta), path(summary) ]
    busco         = BUSCO_BUSCO.out.batch_summary    // channel: [ val(meta), path(summary) ]
    multiqc_files = ch_multiqc_files                 // channel: path(file)
}
