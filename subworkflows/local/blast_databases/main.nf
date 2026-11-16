//
// Build BLAST nucleotide databases (primary, alternate, combined) per assembly.
// --integration_dbs chooses which of those three are built. Mapping, IGV, and
// dotplots consume this channel, so they follow the same set.
//

include { CONCAT_FASTA      } from '../../../modules/local/concat_fasta/main'
include { BLAST_MAKEBLASTDB } from '../../../modules/local/blast/makeblastdb/main'

workflow BLAST_DATABASES {

    take:
    ch_assemblies // channel: [ val(meta), path(primary), path(alternate) ]

    main:
    def mode = params.integration_dbs
    if (!(mode in ['all', 'combined', 'prim_and_alt'])) {
        error("Unknown --integration_dbs '${mode}'. Expected all, combined, or prim_and_alt.")
    }

    // `db_type` distinguishes the databases built for each sample.
    ch_primary = ch_assemblies.map { meta, primary, _alternate -> [meta + [db_type: 'primary'], primary] }
    ch_alt     = ch_assemblies
        .map { meta, _primary, alternate -> [meta + [db_type: 'alternate'], alternate] }
        .filter { _meta, alternate -> alternate }
    ch_haplotypes = ch_primary.mix(ch_alt)

    // Each process is invoked once. prim_and_alt never concatenates; combined
    // never builds the haplotype databases.
    if (mode == 'prim_and_alt') {
        ch_db_input = ch_haplotypes
    } else {
        ch_combined_in = ch_assemblies.map { meta, primary, alternate ->
            [meta + [db_type: 'combined'], [primary, alternate].findAll()]
        }
        CONCAT_FASTA(ch_combined_in)
        ch_db_input = mode == 'combined'
            ? CONCAT_FASTA.out.fasta
            : ch_haplotypes.mix(CONCAT_FASTA.out.fasta)
    }

    BLAST_MAKEBLASTDB(ch_db_input)

    emit:
    db = BLAST_MAKEBLASTDB.out.db // channel: [ val(meta), path(db_files) ]
}
