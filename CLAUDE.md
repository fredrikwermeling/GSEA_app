# GSEA Web Tool — Development Notes

## Project overview
Client-side Gene Set Enrichment Analysis (GSEA) web app.
No backend — all computation runs in the browser via Web Workers.
Hosted on GitHub Pages.

## Tech stack
- Vanilla JavaScript + HTML5 + CSS3
- Plotly.js for figures (CDN)
- PapaParse for CSV/TSV parsing (CDN)
- Web Worker (`worker.js`) for GSEA computation
- No frameworks, no build step

## File structure
- `index.html` — Layout, all CSS (inline), DOM elements
- `app.js` — `GSEAApp` class: UI logic, Plotly rendering, Worker communication
- `worker.js` — GSEA algorithm: enrichment scores, permutations, NES, FDR
- `web_data/` — MSigDB gene set JSON files (H, C2, C5)
- `convert_gmt.py` — Utility to convert GMT → JSON format

## Key conventions
- Design language matches other Wermeling Lab apps (Correlate, Visualize, Green Listed)
- CSS variables defined in `:root` — green-600 (#5a9f4a) is the primary accent
- Cards: green header bar + white body
- Font: Open Sans (body), Roboto Mono (data)

## Gene set JSON format
```json
{ "GENE_SET_NAME": ["GENE1", "GENE2", ...], ... }
```

## Running locally
```bash
cd /Users/fredrikwermeling/Documents/GSEA
python3 -m http.server 8000
# Open http://localhost:8000
```

## Adding new gene set collections
1. Download GMT from MSigDB (https://www.gsea-msigdb.org/)
2. Convert: `python3 convert_gmt.py input.gmt web_data/output.json`
3. Add checkbox in index.html and loading logic in app.js
