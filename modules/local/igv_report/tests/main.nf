include { IGV_REPORT } from '../main'

workflow TEST_IGV_REPORT {
    take:
    reports

    main:
    IGV_REPORT(reports)

    emit:
    report = IGV_REPORT.out.report
    png    = IGV_REPORT.out.png
}
