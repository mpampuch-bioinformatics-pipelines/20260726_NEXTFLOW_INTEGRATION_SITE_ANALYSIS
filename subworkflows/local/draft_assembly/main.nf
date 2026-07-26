//
// Draft genome assembly with hifiasm (primary mode)
//

include { HIFIASM } from '../../../modules/sanger-tol/hifiasm/main'

workflow DRAFT_ASSEMBLY {

    take:
    ch_reads // channel: [ val(meta), path(reads) ]

    main:
    //
    // MODULE: hifiasm assembly. Hi-C, trio and bin/log inputs are left empty as
    //         this pipeline only performs primary/alternate contig assembly.
    //
    HIFIASM(
        ch_reads.map { meta, reads -> [meta, [reads], []] },
        [[:], []],
        [[:], [], []],
        [[:], []],
        [[:], []],
    )

    //
    // Select the primary (p_ctg) and alternate (a_ctg) contig FASTAs, ignoring
    // the unitig/haplotype/balanced graphs also produced by hifiasm.
    //
    ch_assemblies = HIFIASM.out.assembly_fasta
        .map { meta, fastas ->
            def flist     = (fastas instanceof List) ? fastas : [fastas]
            def primary   = flist.find { fa -> fa.name.endsWith('.p_ctg.fa') && !(fa.name =~ /\.(bp|hap[12])\.p_ctg\.fa$/) }
            def alternate = flist.find { fa -> fa.name.endsWith('.a_ctg.fa') }
            [meta, primary, alternate]
        }

    emit:
    assemblies = ch_assemblies // channel: [ val(meta), path(primary), path(alternate) ]
}
