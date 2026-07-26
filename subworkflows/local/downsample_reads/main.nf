//
// Downsample reads to a target coverage based on the reference genome size
//

include { DOWNSAMPLE_SEQKIT } from '../../../modules/local/downsample_seqkit/main'

workflow DOWNSAMPLE_READS {

    take:
    ch_reads_reference // channel: [ val(meta), path(reads), path(reference) ]

    main:
    DOWNSAMPLE_SEQKIT(ch_reads_reference)

    emit:
    reads   = DOWNSAMPLE_SEQKIT.out.reads   // channel: [ val(meta), path(subsampled_reads) ]
    summary = DOWNSAMPLE_SEQKIT.out.summary // channel: [ val(meta), path(tsv) ]
}
