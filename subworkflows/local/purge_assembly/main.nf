//
// Purge haplotypic duplications from the primary assembly with purge_dups
//

include { MINIMAP2_ALIGN as MINIMAP2_READS } from '../../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_SELF  } from '../../../modules/nf-core/minimap2/align/main'
include { PURGEDUPS_PBCSTAT                } from '../../../modules/nf-core/purgedups/pbcstat/main'
include { PURGEDUPS_CALCUTS                } from '../../../modules/nf-core/purgedups/calcuts/main'
include { PURGEDUPS_SPLITFA                } from '../../../modules/nf-core/purgedups/splitfa/main'
include { PURGEDUPS_PURGEDUPS              } from '../../../modules/nf-core/purgedups/purgedups/main'
include { PURGEDUPS_GETSEQS                } from '../../../modules/nf-core/purgedups/getseqs/main'
include { CONCAT_FASTA                     } from '../../../modules/local/concat_fasta/main'

workflow PURGE_ASSEMBLY {

    take:
    ch_assemblies // channel: [ val(meta), path(primary), path(alternate) ]
    ch_reads      // channel: [ val(meta), path(reads) ]

    main:
    //
    // Join reads to their assembly by sample id so alignments are paired correctly.
    // Keep the original alternate for pooling with purge removals later.
    //
    ch_joined = ch_assemblies
        .map { meta, primary, alternate -> [meta.id, meta, primary, alternate] }
        .join(ch_reads.map { meta, reads -> [meta.id, reads] })
        .map { _id, meta, primary, alternate, reads -> [meta, primary, alternate, reads] }

    ch_primary   = ch_joined.map { meta, primary, _alternate, _reads -> [meta, primary] }
    ch_alternate = ch_joined.map { meta, _primary, alternate, _reads -> [meta, alternate] }
    ch_reads_j   = ch_joined.map { meta, _primary, _alternate, reads -> [meta, reads] }

    //
    // MODULE: map reads to the primary assembly and compute base-level coverage
    //
    MINIMAP2_READS(ch_reads_j, ch_primary, false, "", false, false)
    PURGEDUPS_PBCSTAT(MINIMAP2_READS.out.paf)
    PURGEDUPS_CALCUTS(PURGEDUPS_PBCSTAT.out.stat)

    //
    // MODULE: self-align the split primary assembly to detect duplications.
    // An empty reference makes MINIMAP2_ALIGN use the reads as their own target
    // (a self-alignment) while staging the split FASTA only once.
    //
    PURGEDUPS_SPLITFA(ch_primary)
    MINIMAP2_SELF(PURGEDUPS_SPLITFA.out.split_fasta, [[:], []], false, "", false, false)

    //
    // MODULE: identify and remove duplicated sequence
    //
    ch_purge = PURGEDUPS_PBCSTAT.out.basecov
        .join(PURGEDUPS_CALCUTS.out.cutoff)
        .join(MINIMAP2_SELF.out.paf)
    PURGEDUPS_PURGEDUPS(ch_purge)

    ch_getseqs = ch_primary.join(PURGEDUPS_PURGEDUPS.out.bed)
    PURGEDUPS_GETSEQS(ch_getseqs)

    //
    // Pool the original alternate with sequences removed by purge_dups. That
    // pooled set is the alternate assembly for BLAST DBs and hit mapping.
    //
    ch_pool_in = PURGEDUPS_GETSEQS.out.haplotigs
        .join(ch_alternate)
        .map { meta, haplotigs, alternate ->
            [meta + [db_type: 'alternate'], [alternate, haplotigs].findAll()]
        }
    CONCAT_FASTA(ch_pool_in)

    // Strip db_type so downstream joins/stats see the original sample meta.
    ch_pooled_alternate = CONCAT_FASTA.out.fasta.map { meta, fasta ->
        def clean = meta.clone()
        clean.remove('db_type')
        [clean, fasta]
    }

    ch_purged = PURGEDUPS_GETSEQS.out.purged
        .join(ch_pooled_alternate)
        .map { meta, purged, pooled -> [meta, purged, pooled] }

    // One [meta, file] per purge_dups intermediate for sample-centric publishing.
    ch_purgedups_files = PURGEDUPS_PBCSTAT.out.stat
        .mix(PURGEDUPS_PBCSTAT.out.basecov)
        .mix(PURGEDUPS_CALCUTS.out.cutoff)
        .mix(PURGEDUPS_CALCUTS.out.log)
        .mix(PURGEDUPS_SPLITFA.out.split_fasta)
        .mix(PURGEDUPS_PURGEDUPS.out.bed)
        .mix(PURGEDUPS_PURGEDUPS.out.log)

    emit:
    assemblies      = ch_purged           // channel: [ val(meta), path(purged), path(pooled_alternate) ]
    purgedups_files = ch_purgedups_files  // channel: [ val(meta), path(file) ]
}
