//
// Extract BLAST hit subject sequences (with flanks) from the assembly databases,
// and collapse hits on the same contig into one locus window.
//

include { BLAST_EXTRACT_HITS  } from '../../../modules/local/blast/extract_hits/main'
include { BLAST_PREDICT_LOCI } from '../../../modules/local/blast/predict_loci/main'

workflow EXTRACT_SEQUENCES {

    take:
    ch_results // channel: [ val(meta), path(results) ]   (meta carries db_type)
    ch_db      // channel: [ val(meta), path(db_files) ]  (meta carries db_type)

    main:
    // Pair each result file with the database it was searched against (same meta).
    ch_extract_input = ch_db.join(ch_results)
    BLAST_EXTRACT_HITS(ch_extract_input)
    BLAST_PREDICT_LOCI(ch_extract_input)

    emit:
    fasta        = BLAST_EXTRACT_HITS.out.fasta    // channel: [ val(meta), path(extracted_fasta) ]
    all_hits     = BLAST_EXTRACT_HITS.out.hits     // channel: [ val(meta), path(per_hit_fastas) ]
    auto_loci    = BLAST_PREDICT_LOCI.out.loci     // channel: [ val(meta), path(per_locus_fastas) ]
    loci_summary = BLAST_PREDICT_LOCI.out.summary  // channel: [ val(meta), path(loci_tsv) ]
}
