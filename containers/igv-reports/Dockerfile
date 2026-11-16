# igv-reports 1.17.0, samtools, and Playwright's Chromium.
# create_report writes the HTML. Chromium snapshots #igvDiv to PNG because
# igv-reports has no image export (https://github.com/igvteam/igv-reports/issues/78).
# Playwright's image already contains the browser binaries for this tag.
FROM mcr.microsoft.com/playwright/python:v1.49.1-jammy

USER root
ENV DEBIAN_FRONTEND=noninteractive
ENV PIP_DISABLE_PIP_VERSION_CHECK=1

# The base image ships Chromium under /ms-playwright, but a later pip
# install does not keep the matching Playwright Python package. Pin that
# package and do not run `playwright install` (it would download browsers
# we already have).
ENV PLAYWRIGHT_BROWSERS_PATH=/ms-playwright

RUN apt-get update \
    && apt-get install -y --no-install-recommends samtools \
    && pip install --no-cache-dir 'playwright==1.49.1' 'igv-reports==1.17.0' \
    && python -c 'import playwright, igv_reports' \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /data
ENTRYPOINT []
CMD ["/bin/bash"]
