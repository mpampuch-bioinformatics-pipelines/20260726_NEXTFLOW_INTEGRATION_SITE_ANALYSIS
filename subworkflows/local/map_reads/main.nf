//
// Map reads back to the extracted integration-site sequences with minimap2
//

include { MINIMAP2_ALIGN as MINIMAP2_MAP } from '../../../modules/nf-core/minimap2/align/main'

workflow MAP_READS {

    take:
    ch_targets // channel: [ val(meta), path(extracted_fasta) ] (meta carries db_type)
    ch_reads   // channel: [ val(id), path(reads) ]             (per-sample reads)

    main:
    //
    // Pair each extracted-sequence target with its sample's reads.
    //
    ch_map = ch_targets
        .map { meta, target -> [meta.id, meta, target] }
        .combine(ch_reads, by: 0)
        .multiMap { _id, meta, target, reads ->
            reads:  [meta, reads]
            target: [meta, target]
        }

    MINIMAP2_MAP(ch_map.reads, ch_map.target, true, 'bai', false, false)

    emit:
    bam   = MINIMAP2_MAP.out.bam   // channel: [ val(meta), path(bam) ]
    index = MINIMAP2_MAP.out.index // channel: [ val(meta), path(bai) ]
}
