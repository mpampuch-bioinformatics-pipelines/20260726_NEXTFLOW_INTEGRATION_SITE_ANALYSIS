//
// Assembly QC statistics: asmstats, gfastats and BUSCO
//

include { ASMSTATS              } from '../../../modules/sanger-tol/asmstats/main'
include { GFASTATS              } from '../../../modules/nf-core/gfastats/main'
include { BUSCO_BUSCO           } from '../../../modules/nf-core/busco/busco/main'
include { RENAME_BUSCO_SUMMARY  } from '../../../modules/local/rename_busco_summary/main'

workflow GENOME_STATISTICS {

    take:
    ch_assemblies  // channel: [ val(meta), path(primary), path(alternate) ]
    val_stage      // string: assembly stage label (e.g. 'raw' or 'purged')
    ch_lineage_db  // channel: [ val(lineage), path(busco_downloads) ]

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
    // MODULE: gene-completeness assessment with BUSCO.
    // A named lineage reuses the shared busco/download directory. Auto-lineage
    // still picks a dataset from the assembly, so those tasks download their own.
    //
    ch_busco_keyed = ch_split.map { meta, fasta ->
        def lineage = meta.busco_lineage ?: params.busco_lineage ?: 'auto_euk'
        [lineage, meta + [busco_lineage: lineage], fasta]
    }
    // Auto-lineage has no shared download. An empty path tells BUSCO to pick one.
    ch_busco_auto = ch_busco_keyed
        .filter { lineage, _meta, _fasta -> lineage in ['auto', 'auto_prok', 'auto_euk'] }
        .map { _lineage, meta, fasta -> [meta, fasta, []] }
    ch_busco_named = ch_busco_keyed
        .filter { lineage, _meta, _fasta -> !(lineage in ['auto', 'auto_prok', 'auto_euk']) }
        .combine(ch_lineage_db, by: 0)
        .map { _lineage, meta, fasta, db -> [meta, fasta, db] }
    // multiMap keeps each assembly aligned with its own lineage directory.
    ch_busco = ch_busco_auto.mix(ch_busco_named).multiMap { meta, fasta, db ->
        fasta: [meta, fasta]
        lineages: db
    }

    BUSCO_BUSCO(
        ch_busco.fasta,
        'genome',
        params.busco_lineage ?: 'auto_euk',
        ch_busco.lineages,
        [],
        true,
    )
    // short_summary*.txt is shared across haplotypes until the sample id is inserted.
    RENAME_BUSCO_SUMMARY(BUSCO_BUSCO.out.short_summaries_txt)
    ch_multiqc_files = ch_multiqc_files.mix(RENAME_BUSCO_SUMMARY.out.summary.map { _meta, summary -> summary })

    emit:
    asmstats      = ASMSTATS.out.stats               // channel: [ val(meta), path(stats) ]
    gfastats      = GFASTATS.out.assembly_summary    // channel: [ val(meta), path(summary) ]
    busco         = BUSCO_BUSCO.out.batch_summary    // channel: [ val(meta), path(summary) ]
    multiqc_files = ch_multiqc_files                 // channel: path(file)
}
