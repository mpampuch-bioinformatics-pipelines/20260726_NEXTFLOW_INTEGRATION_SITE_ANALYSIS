//
// Lift reference-genome genes onto each extracted hit.
// Runs only for samples that carry an annotation file.
//

include { LIFT_REFERENCE_ANNOTATIONS as LIFT_ONE_HIT } from '../../../modules/local/lift_reference_annotations/main'

workflow LIFT_REFERENCE_ANNOTATIONS {
    take:
    ch_hits       // channel: [ val(meta), path(hit_fastas) ]  (hit_fastas may be a list; meta carries db_type)
    ch_reference  // channel: [ val(meta), path(reference) ]
    ch_annotation // channel: [ val(id), path(annotation) ]

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

    // Inner join: a sample with no annotation is not lifted.
    ch_input = ch_per_hit
        .combine(ch_ref, by: 0)
        .combine(ch_annotation, by: 0)
        .map { _id, meta, hit, reference, annotation ->
            [meta, hit, reference, annotation]
        }

    LIFT_ONE_HIT(ch_input)

    emit:
    lifted     = LIFT_ONE_HIT.out.lifted_gff3.mix(LIFT_ONE_HIT.out.lifted_gtf)
    validation = LIFT_ONE_HIT.out.validation
    unmapped   = LIFT_ONE_HIT.out.unmapped
}
