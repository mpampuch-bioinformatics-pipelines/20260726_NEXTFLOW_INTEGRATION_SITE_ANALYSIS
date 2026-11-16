//
// Measure the four manual junction checks on the read-inspector BAM.
// The table does not choose trust or revisit.
//

include { JUNCTION_QC as SCORE_JUNCTIONS } from '../../../modules/local/junction_qc/main'

workflow JUNCTION_QC {
    take:
    ch_inspector // channel: [ val(meta), path(dir) ]   (read-inspector directory; meta carries db_type)
    ch_loci      // channel: [ val(meta), path(tsv) ]   (loci.tsv or hit_intervals.tsv)

    main:
    // Join on sample and assembly database. The meta maps are not identical:
    // the inspector directory and the locus table are separate processes.
    ch_input = ch_inspector
        .map { meta, inspector ->
            ["${meta.id}|${meta.db_type}", meta, inspector]
        }
        .combine(
            ch_loci.map { meta, loci ->
                ["${meta.id}|${meta.db_type}", loci]
            },
            by: 0
        )
        .map { _key, meta, inspector, loci ->
            [meta, inspector, loci]
        }

    SCORE_JUNCTIONS(ch_input)

    emit:
    metrics = SCORE_JUNCTIONS.out.metrics // channel: [ val(meta), path(tsv) ]
}
