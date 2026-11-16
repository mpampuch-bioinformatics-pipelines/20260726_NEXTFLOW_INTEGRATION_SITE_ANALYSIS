//
// Extract plasmid features and annotate extracted integration-site hit windows.
//

include { PLASMID_EXTRACT_FEATURES } from '../../../modules/local/plasmid/extract_features/main'
include { PLASMID_BLAST_FEATURES   } from '../../../modules/local/plasmid/blast_features/main'
include { PLASMID_ANNOTATE_HIT     } from '../../../modules/local/plasmid/annotate_hit/main'

// Same plasmid file is extracted once. The id is the file name plus a hash of
// its path so two different plasmids with the same name do not collide.
def sharedMeta(String key, String label) {
    def hash = java.security.MessageDigest.getInstance('MD5')
        .digest(key.getBytes('UTF-8'))
        .encodeHex()
        .toString()
        .substring(0, 8)
    def base = key.contains('|') ? label : new File(key).name
    def safe = base.replaceAll(/[^A-Za-z0-9._-]/, '_')
    [id: "${safe}_${hash}", key: key]
}

workflow PLASMID_ANNOTATION {

    take:
    ch_plasmid // channel: [ val(meta), path(plasmid) ]  (only samples with a plasmid)
    ch_hits    // channel: [ val(meta), path(hit_fastas) ]  (meta carries db_type)

    main:
    //
    // One extraction per plasmid path, then copy that result onto each sample.
    //
    ch_keyed = ch_plasmid.map { meta, plasmid ->
        def plasmid_file = file(plasmid, checkIfExists: true)
        [plasmid_file.toAbsolutePath().toString(), meta, plasmid_file]
    }
    PLASMID_EXTRACT_FEATURES(
        ch_keyed
            .map { key, _meta, plasmid_file -> [sharedMeta(key, 'plasmid'), plasmid_file] }
            .unique { shared, _plasmid -> shared.key }
    )
    ch_sample_by_key = ch_keyed.map { key, meta, _plasmid -> [key, meta] }

    ch_features_by_id = PLASMID_EXTRACT_FEATURES.out.features_fasta
        .map { shared, fasta -> [shared.key, fasta] }
        .combine(ch_sample_by_key, by: 0)
        .map { _key, fasta, meta -> [meta.id, fasta] }

    ch_plasmid_gb = PLASMID_EXTRACT_FEATURES.out.gb
        .map { shared, gb -> [shared.key, gb] }
        .combine(ch_sample_by_key, by: 0)
        .map { _key, gb, meta -> [meta, gb] }
    ch_plasmid_gff3 = PLASMID_EXTRACT_FEATURES.out.gff3
        .map { shared, gff3 -> [shared.key, gff3] }
        .combine(ch_sample_by_key, by: 0)
        .map { _key, gff3, meta -> [meta, gff3] }
    ch_plasmid_features = PLASMID_EXTRACT_FEATURES.out.features_fasta
        .map { shared, fasta -> [shared.key, fasta] }
        .combine(ch_sample_by_key, by: 0)
        .map { _key, fasta, meta -> [meta, fasta] }
    // Qualifiers (product, note, codon_start) are not used by BLAST.
    ch_features_meta_by_id = PLASMID_EXTRACT_FEATURES.out.features_meta
        .map { shared, features_meta -> [shared.key, features_meta] }
        .combine(ch_sample_by_key, by: 0)
        .map { _key, features_meta, meta -> [meta.id, features_meta] }

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
    ch_annotate = PLASMID_BLAST_FEATURES.out.blast
        .map { meta, hit, features, blast -> [meta.id, meta, hit, features, blast] }
        .combine(ch_features_meta_by_id, by: 0)
        .map { _id, meta, hit, features, blast, features_meta ->
            [meta, hit, features, blast, features_meta]
        }
    PLASMID_ANNOTATE_HIT(ch_annotate)

    emit:
    plasmid_gb             = ch_plasmid_gb        // channel: [ val(meta), path(gb) ]
    plasmid_gff3           = ch_plasmid_gff3      // channel: [ val(meta), path(gff3) ]
    plasmid_features_fasta = ch_plasmid_features  // channel: [ val(meta), path(fasta) ]
    feature_blast          = PLASMID_BLAST_FEATURES.out.blast.map { meta, _hit, _feat, blast -> [meta, blast] }
    hit_gff3               = PLASMID_ANNOTATE_HIT.out.gff3               // channel: [ val(meta), path(gff3) ]
    hit_gb                 = PLASMID_ANNOTATE_HIT.out.gb                 // channel: [ val(meta), path(gb) ]
}
