//
// Extract plasmid features and annotate extracted integration-site hit windows.
//

include { PLASMID_EXTRACT_FEATURES } from '../../../modules/local/plasmid/extract_features/main'
include { PLASMID_BLAST_FEATURES   } from '../../../modules/local/plasmid/blast_features/main'
include { PLASMID_ANNOTATE_HIT     } from '../../../modules/local/plasmid/annotate_hit/main'

workflow PLASMID_ANNOTATION {

    take:
    ch_plasmid // channel: [ val(meta), path(plasmid) ]  (only samples with a plasmid)
    ch_hits    // channel: [ val(meta), path(hit_fastas) ]  (meta carries db_type)

    main:
    //
    // Normalize SnapGene/GenBank → feature GFF3 + feature FASTA (+ GenBank copy).
    //
    PLASMID_EXTRACT_FEATURES(ch_plasmid)

    //
    // Join per-sample feature FASTA to extracted hit windows, then one job per hit.
    //
    ch_features_by_id = PLASMID_EXTRACT_FEATURES.out.features_fasta.map { meta, fasta ->
        [meta.id, fasta]
    }

    ch_per_hit = ch_hits
        .map { meta, hits ->
            def hit_list = hits instanceof List ? hits : [hits]
            [meta, hit_list]
        }
        .transpose()
        .map { meta, hit ->
            [meta.id, meta + [hit_id: hit.baseName], hit]
        }
        .combine(ch_features_by_id, by: 0)
        .map { _id, meta, hit, features ->
            [meta, hit, features]
        }

    PLASMID_BLAST_FEATURES(ch_per_hit)
    PLASMID_ANNOTATE_HIT(PLASMID_BLAST_FEATURES.out.blast)

    emit:
    plasmid_gb             = PLASMID_EXTRACT_FEATURES.out.gb             // channel: [ val(meta), path(gb) ]
    plasmid_gff3           = PLASMID_EXTRACT_FEATURES.out.gff3           // channel: [ val(meta), path(gff3) ]
    plasmid_features_fasta = PLASMID_EXTRACT_FEATURES.out.features_fasta // channel: [ val(meta), path(fasta) ]
    feature_blast          = PLASMID_BLAST_FEATURES.out.blast.map { meta, _hit, _feat, blast -> [meta, blast] }
    hit_gff3               = PLASMID_ANNOTATE_HIT.out.gff3               // channel: [ val(meta), path(gff3) ]
    hit_gb                 = PLASMID_ANNOTATE_HIT.out.gb                 // channel: [ val(meta), path(gb) ]
}
