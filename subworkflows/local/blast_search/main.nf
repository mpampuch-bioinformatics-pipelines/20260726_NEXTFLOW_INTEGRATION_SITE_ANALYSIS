//
// BLAST the combined construct queries against each assembly database
//

include { BLAST_BLASTN } from '../../../modules/local/blast/blastn/main'

workflow BLAST_SEARCH {

    take:
    ch_db    // channel: [ val(meta), path(db_files) ]  (meta carries db_type)
    ch_query // channel: [ val(id), path(query_fasta) ] (per-sample merged query)

    main:
    //
    // Attach the per-sample query FASTA to each of that sample's databases.
    //
    ch_blastn_input = ch_db
        .map { meta, db -> [meta.id, meta, db] }
        .combine(ch_query, by: 0)
        .map { _id, meta, db, query -> [meta, db, query] }

    BLAST_BLASTN(ch_blastn_input)

    emit:
    results = BLAST_BLASTN.out.results // channel: [ val(meta), path(results) ]
}
