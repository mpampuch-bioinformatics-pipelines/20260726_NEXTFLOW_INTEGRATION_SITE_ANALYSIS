#!/usr/bin/env python3
"""Save PNG snapshots of an igv-reports HTML file.

create_report does not write images. The HTML page assigns the browser to
`igvBrowser` and draws it in `#igvDiv`. One PNG per variant-table row, which
is one locus for this pipeline. A single row is named `<prefix>_igv_mqc.png`
so MultiQC embeds it. Further rows are `<prefix>_igv_<n>_mqc.png`.
"""

import asyncio
import sys
from pathlib import Path

from playwright.async_api import async_playwright


async def snapshot(html_path, prefix):
    html = Path(html_path).resolve()
    if not html.is_file():
        raise SystemExit(f"IGV report not found: {html}")

    async with async_playwright() as playwright:
        browser = await playwright.chromium.launch(
            args=["--no-sandbox", "--disable-dev-shm-usage"]
        )
        page = await browser.new_page(viewport={"width": 1400, "height": 900})
        await page.goto(html.as_uri(), wait_until="domcontentloaded")
        # `let igvBrowser` is global to the page script but not window.igvBrowser.
        await page.wait_for_function(
            "() => typeof igvBrowser !== 'undefined' && igvBrowser",
            timeout=120000,
        )
        await page.wait_for_selector("#igvDiv canvas", timeout=120000)
        await page.wait_for_timeout(1500)

        rows = await page.query_selector_all("#variant_table tbody tr")
        count = len(rows) if rows else 1
        if not rows:
            rows = [None]

        written = []
        for index, row in enumerate(rows, start=1):
            if row is not None and count > 1:
                await row.click()
                await page.wait_for_selector("#igvDiv canvas", timeout=120000)
                await page.wait_for_timeout(1500)
            name = (
                f"{prefix}_igv_mqc.png"
                if count == 1
                else f"{prefix}_igv_{index}_mqc.png"
            )
            panel = await page.query_selector("#igvDiv")
            if panel is None:
                raise SystemExit("IGV report has no #igvDiv panel to snapshot")
            await panel.screenshot(path=name)
            written.append(name)

        await browser.close()

    if not written:
        raise SystemExit("IGV snapshot wrote no PNG")


def main():
    if len(sys.argv) != 3:
        raise SystemExit(f"usage: {sys.argv[0]} REPORT.html PREFIX")
    asyncio.run(snapshot(sys.argv[1], sys.argv[2]))


if __name__ == "__main__":
    main()
