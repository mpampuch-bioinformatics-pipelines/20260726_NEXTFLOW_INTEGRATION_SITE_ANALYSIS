//
// Map the original HiFi reads to every selected hit of one assembly and
// write the IGV / JBrowse directory. This is not the per-hit HTML report.
//

include { READ_INSPECTOR as INSPECT_READS } from '../../../modules/local/read_inspector/main'

workflow READ_INSPECTOR {
    take:
    ch_hits  // channel: [ val(meta), path(hit_fastas) ]  (one element per sample and db_type; hits may be a list)
    ch_reads // channel: [ val(id), path(reads) ]         (original sample FASTQ, not the downsampled reads)

    main:
    ch_input = ch_hits
        .map { meta, hits ->
            [meta.id, meta, hits]
        }
        .combine(ch_reads, by: 0)
        .map { _id, meta, hits, reads ->
            [meta, reads, hits]
        }

    INSPECT_READS(ch_input)

    emit:
    inspector = INSPECT_READS.out.inspector // channel: [ val(meta), path(dir) ]
}
