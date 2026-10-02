#!/usr/bin/env bash -C -e -u -o pipefail
mkdir multiqc_data
mkdir multiqc_plots
touch multiqc_report.html

cat <<-END_VERSIONS > versions.yml
"NFCORE_INTEGRATIONSITEANALYSIS:INTEGRATION_SITE_ANALYSIS:MULTIQC":
    multiqc: $( multiqc --version | sed -e "s/multiqc, version //g" )
END_VERSIONS
