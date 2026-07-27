//
// Read QC: kmer profiling (FastK + GenomeScope2) and read summary stats (seqkit)
//

include { SEQKIT_STATS } from '../../../modules/local/seqkit_stats/main'
include { FASTK_FASTK   } from '../../../modules/nf-core/fastk/fastk/main'
include { FASTK_HISTEX  } from '../../../modules/nf-core/fastk/histex/main'
include { GENOMESCOPE2  } from '../../../modules/nf-core/genomescope2/main'

workflow READ_QC {

    take:
    ch_reads // channel: [ val(meta), path(reads) ]

    main:
    ch_multiqc_files = channel.empty()

    //
    // MODULE: Summary statistics of the input reads
    //
    SEQKIT_STATS(ch_reads)
    ch_multiqc_files = ch_multiqc_files.mix(SEQKIT_STATS.out.stats.map { _meta, stats -> stats })

    //
    // MODULE: Build a FastK kmer database and estimate genome characteristics
    //
    FASTK_FASTK(ch_reads)
    FASTK_HISTEX(FASTK_FASTK.out.hist)
    GENOMESCOPE2(FASTK_HISTEX.out.hist)
    ch_multiqc_files = ch_multiqc_files.mix(GENOMESCOPE2.out.summary.map { _meta, summary -> summary })

    // One [meta, file] per GenomeScope artifact for sample-centric publishing.
    ch_genomescope_plots = GENOMESCOPE2.out.linear_plot_png
        .mix(GENOMESCOPE2.out.transformed_linear_plot_png)
        .mix(GENOMESCOPE2.out.log_plot_png)
        .mix(GENOMESCOPE2.out.transformed_log_plot_png)
        .mix(GENOMESCOPE2.out.model)
        .mix(GENOMESCOPE2.out.json_report)

    emit:
    read_stats        = SEQKIT_STATS.out.stats  // channel: [ val(meta), path(tsv) ]
    genomescope       = GENOMESCOPE2.out.summary // channel: [ val(meta), path(summary) ]
    genomescope_plots = ch_genomescope_plots     // channel: [ val(meta), path(file) ]
    multiqc_files     = ch_multiqc_files         // channel: path(file)
}
