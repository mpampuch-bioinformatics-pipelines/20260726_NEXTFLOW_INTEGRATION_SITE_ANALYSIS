include { RENAME_BUSCO_SUMMARY } from '../main'

workflow TEST_RENAME_BUSCO_SUMMARY {
    take:
    summaries

    main:
    RENAME_BUSCO_SUMMARY(summaries)

    emit:
    summary = RENAME_BUSCO_SUMMARY.out.summary
}
