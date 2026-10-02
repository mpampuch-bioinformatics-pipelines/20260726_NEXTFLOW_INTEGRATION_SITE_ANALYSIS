//
// Downsample reads to a target coverage. Genome size is a shared seqkit/stats TSV.
//

include { DOWNSAMPLE_SEQKIT } from '../../../modules/local/downsample_seqkit/main'

workflow DOWNSAMPLE_READS {

    take:
    ch_reads_size // channel: [ val(meta), path(reads), path(genome_size_tsv) ]

    main:
    DOWNSAMPLE_SEQKIT(ch_reads_size)

    emit:
    reads   = DOWNSAMPLE_SEQKIT.out.reads   // channel: [ val(meta), path(subsampled_reads) ]
    summary = DOWNSAMPLE_SEQKIT.out.summary // channel: [ val(meta), path(tsv) ]
}
