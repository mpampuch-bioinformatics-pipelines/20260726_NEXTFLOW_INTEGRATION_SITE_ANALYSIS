//
// Build BLAST nucleotide databases (primary, alternate, combined) per assembly
//

include { CONCAT_FASTA      } from '../../../modules/local/concat_fasta/main'
include { BLAST_MAKEBLASTDB } from '../../../modules/local/blast/makeblastdb/main'

workflow BLAST_DATABASES {

    take:
    ch_assemblies // channel: [ val(meta), path(primary), path(alternate) ]

    main:
    // `db_type` distinguishes the three databases built for each sample.
    ch_primary = ch_assemblies.map { meta, primary, _alternate -> [meta + [db_type: 'primary'], primary] }
    ch_alt     = ch_assemblies
        .map { meta, _primary, alternate -> [meta + [db_type: 'alternate'], alternate] }
        .filter { _meta, alternate -> alternate }

    //
    // MODULE: concatenate primary + alternate into a single combined FASTA
    //
    ch_combined_in = ch_assemblies.map { meta, primary, alternate -> [meta + [db_type: 'combined'], [primary, alternate].findAll()] }
    CONCAT_FASTA(ch_combined_in)

    //
    // MODULE: build a BLAST database for each FASTA
    //
    ch_db_input = ch_primary.mix(ch_alt).mix(CONCAT_FASTA.out.fasta)
    BLAST_MAKEBLASTDB(ch_db_input)

    emit:
    db = BLAST_MAKEBLASTDB.out.db // channel: [ val(meta), path(db_files) ]
}
