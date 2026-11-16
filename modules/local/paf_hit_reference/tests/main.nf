include { PAF_HIT_REFERENCE } from '../main'

workflow TEST_PAF_HIT_REFERENCE {
    take:
    alignments

    main:
    PAF_HIT_REFERENCE(alignments)

    emit:
    best_chromosome = PAF_HIT_REFERENCE.out.best_chromosome
    best_fasta      = PAF_HIT_REFERENCE.out.best_fasta
}
