#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PLASMID_ANNOTATION } from '../subworkflows/local/plasmid_annotation/main'

workflow {
    def root = file("${projectDir}/..")
    ch_plasmid = Channel.of([
        [id: 'test'],
        file("${root}/TESTS/TEST_DATA/PLASMIDS/cmc_087.gb", checkIfExists: true)
    ])
    ch_hits = Channel.of([
        [id: 'test', db_type: 'primary'],
        [
            file("${root}/TESTS/TEST_DATA/PLASMIDS/stub_hit.fa", checkIfExists: true)
        ]
    ])

    PLASMID_ANNOTATION(ch_plasmid, ch_hits)

    PLASMID_ANNOTATION.out.hit_gff3.view { meta, gff -> "hit_gff3: ${meta.hit_id} -> ${gff}" }
    PLASMID_ANNOTATION.out.hit_gb.view { meta, gb -> "hit_gb: ${meta.hit_id} -> ${gb}" }
}
