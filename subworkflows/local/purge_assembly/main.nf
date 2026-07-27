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

workflow PURGE_ASSEMBLY {

    take:
    ch_assemblies // channel: [ val(meta), path(primary), path(alternate) ]
    ch_reads      // channel: [ val(meta), path(reads) ]

    main:
    //
    // Join reads to their assembly by sample id so alignments are paired correctly.
    //
    ch_joined = ch_assemblies
        .map { meta, primary, alternate -> [meta.id, meta, primary, alternate] }
        .join(ch_reads.map { meta, reads -> [meta.id, reads] })
        .map { _id, meta, primary, _alternate, reads -> [meta, primary, reads] }

    ch_primary = ch_joined.map { meta, primary, _reads -> [meta, primary] }
    ch_reads_j = ch_joined.map { meta, _primary, reads -> [meta, reads] }

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
    // The purged primary becomes the new primary; retained haplotigs become the
    // new "alternate" assembly for downstream BLAST database construction.
    //
    ch_purged = PURGEDUPS_GETSEQS.out.purged
        .join(PURGEDUPS_GETSEQS.out.haplotigs)
        .map { meta, purged, haplotigs -> [meta, purged, haplotigs] }

    // One [meta, file] per purge_dups intermediate for sample-centric publishing.
    ch_purgedups_files = PURGEDUPS_PBCSTAT.out.stat
        .mix(PURGEDUPS_PBCSTAT.out.basecov)
        .mix(PURGEDUPS_CALCUTS.out.cutoff)
        .mix(PURGEDUPS_CALCUTS.out.log)
        .mix(PURGEDUPS_SPLITFA.out.split_fasta)
        .mix(PURGEDUPS_PURGEDUPS.out.bed)
        .mix(PURGEDUPS_PURGEDUPS.out.log)

    emit:
    assemblies      = ch_purged           // channel: [ val(meta), path(purged), path(haplotigs) ]
    purgedups_files = ch_purgedups_files  // channel: [ val(meta), path(file) ]
}
