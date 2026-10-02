//
// Assembly QC with QUAST: compare assemblies to reference genome
//

include { QUAST } from '../../../modules/nf-core/quast/main'

workflow ASSEMBLY_QC {

    take:
    ch_assemblies // channel: [ val(meta), path(primary), path(alternate) ]
    ch_reference  // channel: [ val(meta), path(reference) ]
    val_stage     // string: assembly stage label (e.g. 'raw' or 'purged')

    main:
    ch_multiqc_files = channel.empty()

    //
    // Split primary and alternate assemblies per haplotype and stamp meta.stage / meta.assembly
    //
    ch_split = ch_assemblies
        .flatMap { meta, primary, alternate ->
            def meta_stage = meta + [stage: val_stage]
            def out = [[meta_stage + [assembly: 'primary'], primary]]
            if (alternate) {
                out << [meta_stage + [assembly: 'alternate'], alternate]
            }
            return out
        }
        .filter { _meta, assembly -> assembly }

    //
    // Pair the sample reference with every haplotype. join() would consume the
    // one reference on the first haplotype and drop primary or alternate.
    //
    ch_ref_by_id = ch_reference.map { meta_ref, ref ->
        def id = (meta_ref instanceof Map) ? meta_ref.id : meta_ref
        [id, ref]
    }

    ch_quast_input = ch_split
        .map { meta, assembly -> [meta.id, meta, assembly] }
        .combine(ch_ref_by_id, by: 0)
        .map { _id, meta, assembly, reference ->
            def ref_meta = [id: meta.id]
            def ref_path = reference ? (reference instanceof List ? reference : [reference]) : []
            [ meta, assembly, ref_meta, ref_path, [id: 'gff'], [] ]
        }

    //
    // MODULE: QUAST
    //
    QUAST (
        ch_quast_input.map { meta, assembly, _ref_meta, _ref_path, _gff_meta, _gff_path -> [meta, assembly] },
        ch_quast_input.map { _meta, _assembly, ref_meta, ref_path, _gff_meta, _gff_path -> [ref_meta, ref_path] },
        ch_quast_input.map { _meta, _assembly, _ref_meta, _ref_path, gff_meta, gff_path -> [gff_meta, gff_path] }
    )

    ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.tsv.map { _meta, tsv -> tsv })

    emit:
    results       = QUAST.out.results        // channel: [ val(meta), path(dir) ]
    tsv           = QUAST.out.tsv            // channel: [ val(meta), path(tsv) ]
    multiqc_files = ch_multiqc_files          // channel: path(file)
}
