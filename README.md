# GSEA Web Tool

**Live demo: [https://fredrikwermeling.github.io/GSEA_app/](https://fredrikwermeling.github.io/GSEA_app/)**

Part of the Wermeling Lab tool family with [Green Listed](https://greenlisted.cmm.se) and [Correlate](https://correlate.cmm.se); the three apps share one palette and header layout.

A client-side Gene Set Enrichment Analysis (GSEA) web application. All computation runs in the browser — no backend required.

## Features

- **Preranked GSEA** — weighted enrichment score (fgsea-style) with gene-label permutation testing
- **MSigDB gene sets** — Hallmark (H), Curated (C2), GO (C5) collections included, plus custom GMT upload
- **Interactive figures** — bubble/dot plot, running enrichment score curve with rug plot, ranked gene list
- **Publication-quality export** — SVG and PNG at configurable resolution
- **Results table** — sortable, filterable, downloadable as CSV
- **Auto-generated methods section** — copy-paste ready for manuscripts
- **Web Worker computation** — UI stays responsive during analysis

## Usage

1. Upload a CSV or TSV file with gene-level statistics (e.g. DESeq2 or MAGeCK output)
2. Select the gene name column and ranking metric column
3. Choose gene set collections and analysis parameters
4. Click **Run GSEA**

Or click **Load Example Data** to try it with simulated data.

## Running locally

```bash
git clone https://github.com/fredrikwermeling/GSEA_app.git
cd GSEA_app
python3 -m http.server 8000
# Open http://localhost:8000
```

## Gene set collections

| Collection | Sets | Description |
|-----------|------|-------------|
| H: Hallmark | 50 | Well-defined biological states and processes |
| C2: Curated | 7,233 | KEGG, Reactome, WikiPathways and more |
| C5: GO | 16,008 | Gene Ontology BP, MF, CC |

Gene sets from [MSigDB](https://www.gsea-msigdb.org/) v2023.2. Custom GMT files can also be uploaded.

## Logo and design assets

- `logo-word.svg` — outline-traced "Enrich" wordmark (Open Sans 700, fill `#83AA3D`, the same face and green as the Green Listed wordmark). Regenerate with fontTools from `greenlistedv2/fonts/open-sans-latin.woff2` if the name changes.
- `logo-mark.svg` — running-enrichment-score mark shown to the right of the wordmark.
- `favicon.svg` / `favicon.png` — bolder variant of the mark.
- Colour tokens live in `:root` in `index.html` and are copied verbatim from Green Listed / Correlate (accent `--green-600: #6ba544`).

## Adding custom gene set collections

1. Download a GMT file from [MSigDB](https://www.gsea-msigdb.org/)
2. Convert to JSON: `python3 convert_gmt.py input.gmt web_data/output.json`
3. Add a checkbox and loading logic in the app

## Tech stack

- Vanilla JavaScript — no frameworks, no build step
- [Plotly.js](https://plotly.com/javascript/) for figures
- [PapaParse](https://www.papaparse.com/) for CSV/TSV parsing
- Web Workers for GSEA computation

## References

- Subramanian et al., *Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles.* PNAS, 2005.
- Korotkevich et al., *Fast gene set enrichment analysis.* bioRxiv, 2021.
- Liberzon et al., *The Molecular Signatures Database (MSigDB) hallmark gene set collection.* Cell Systems, 2015.

## License

MIT

---

*Developed by [Wermeling Lab](https://wermelinglab.com), Karolinska Institutet*
