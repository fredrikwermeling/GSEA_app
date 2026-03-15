#!/usr/bin/env python3
"""
Prepare DepMap example datasets for Enrich web app.

Downloads expression and CRISPR data from DepMap for 5 cell lines,
computes per-cell-line ranking metrics, and saves as JSON files.

Expression: z-score vs median across all DepMap cell lines
CRISPR: raw Chronos gene effect scores (negative = dependency)
"""

import json
import os
import sys

try:
    import pandas as pd
    import numpy as np
except ImportError:
    print("Install dependencies: pip install pandas numpy")
    sys.exit(1)

# Cell lines to process
CELL_LINES = {
    'A375':  'ACH-000219',  # Skin cancer (Melanoma)
    'A549':  'ACH-000681',  # Lung cancer (Adenocarcinoma)
    'HT29':  'ACH-000552',  # Colon cancer (Colorectal)
    'Raji':  'ACH-000007',  # Blood cancer (B-cell Lymphoma)
    'U251':  'ACH-000232',  # Brain cancer (Glioblastoma)
}

OUTPUT_DIR = 'web_data'

# DepMap file URLs (25Q3 release)
EXPRESSION_URL = 'https://depmap.org/portal/download/api/download?file_name=downloads-by-canonical-id%2Fpublic-25q3-b56c.80%2FOmicsExpressionTPMLogp1HumanProteinCodingGenes.csv&dl_name=OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv&bucket=depmap-external-downloads'
CRISPR_URL = 'https://depmap.org/portal/download/api/download?file_name=downloads-by-canonical-id%2F25q3-public-6202.1%2FCRISPRGeneEffect.csv&dl_name=CRISPRGeneEffect.csv&bucket=depmap-external-downloads'


def clean_gene_name(col_name):
    """Extract gene symbol from DepMap column format: 'GENE (12345)'"""
    if '(' in col_name:
        return col_name.split('(')[0].strip()
    return col_name.strip()


def download_if_needed(url, local_path):
    """Download file if not cached locally."""
    if os.path.exists(local_path):
        print(f"  Using cached: {local_path}")
        return
    print(f"  Downloading: {url}")
    print(f"  (This may take a few minutes for large files...)")
    import urllib.request
    urllib.request.urlretrieve(url, local_path)
    print(f"  Saved to: {local_path}")


def process_expression(df, cell_lines):
    """Compute z-scores for each cell line vs median across all cell lines."""
    # Compute median and MAD across all cell lines
    median = df.median(axis=0)
    mad = (df - median).abs().median(axis=0)
    # Use 5th percentile of non-zero MADs as floor (avoids extreme z-scores
    # from near-zero variability genes while preserving true signal)
    nonzero_mads = mad[mad > 0]
    min_mad = nonzero_mads.quantile(0.05) if len(nonzero_mads) > 0 else 0.1
    mad = mad.clip(lower=min_mad)

    results = {}
    for name, ach_id in cell_lines.items():
        if ach_id not in df.index:
            print(f"  WARNING: {name} ({ach_id}) not found in expression data!")
            continue
        row = df.loc[ach_id]
        zscores = (row - median) / mad
        # Build gene list, sorted by z-score descending
        genes = []
        for col in df.columns:
            gene = clean_gene_name(col)
            z = round(float(zscores[col]), 4)
            if not np.isnan(z):
                genes.append({'Gene': gene, 'Expression_zscore': z})
        genes.sort(key=lambda x: x['Expression_zscore'], reverse=True)
        results[name] = genes
        print(f"  {name}: {len(genes)} genes, range [{genes[-1]['Expression_zscore']:.2f}, {genes[0]['Expression_zscore']:.2f}]")
    return results


def process_crispr(df, cell_lines):
    """Extract raw Chronos gene effect scores for each cell line."""
    results = {}
    for name, ach_id in cell_lines.items():
        if ach_id not in df.index:
            print(f"  WARNING: {name} ({ach_id}) not found in CRISPR data!")
            continue
        row = df.loc[ach_id]
        genes = []
        for col in df.columns:
            gene = clean_gene_name(col)
            score = round(float(row[col]), 4) if not np.isnan(row[col]) else None
            if score is not None:
                genes.append({'Gene': gene, 'Chronos_score': score})
        genes.sort(key=lambda x: x['Chronos_score'], reverse=True)
        results[name] = genes
        print(f"  {name}: {len(genes)} genes, range [{genes[-1]['Chronos_score']:.2f}, {genes[0]['Chronos_score']:.2f}]")
    return results


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    cache_dir = os.path.join(OUTPUT_DIR, '.cache')
    os.makedirs(cache_dir, exist_ok=True)

    expr_cache = os.path.join(cache_dir, 'OmicsExpressionProteinCodingGenesTPMLogp1.csv')
    crispr_cache = os.path.join(cache_dir, 'CRISPRGeneEffect.csv')

    # Download data
    print("Step 1: Downloading DepMap data...")
    download_if_needed(EXPRESSION_URL, expr_cache)
    download_if_needed(CRISPR_URL, crispr_cache)

    # Process expression — file has metadata columns before gene columns
    print("\nStep 2: Processing expression data...")
    expr_raw = pd.read_csv(expr_cache, index_col=0)
    # Use ModelID as index, drop non-gene metadata columns
    meta_cols = ['SequencingID', 'ModelID', 'IsDefaultEntryForModel',
                 'ModelConditionID', 'IsDefaultEntryForMC']
    existing_meta = [c for c in meta_cols if c in expr_raw.columns]
    if 'ModelID' in expr_raw.columns:
        expr_df = expr_raw.set_index('ModelID').drop(columns=[c for c in existing_meta if c != 'ModelID'], errors='ignore')
        # Keep only numeric columns
        expr_df = expr_df.select_dtypes(include=[np.number])
    else:
        expr_df = expr_raw.select_dtypes(include=[np.number])
    print(f"  Expression matrix: {expr_df.shape[0]} cell lines x {expr_df.shape[1]} genes")
    expr_results = process_expression(expr_df, CELL_LINES)

    # Process CRISPR
    print("\nStep 3: Processing CRISPR data...")
    crispr_df = pd.read_csv(crispr_cache, index_col=0)
    crispr_results = process_crispr(crispr_df, CELL_LINES)

    # Save JSON files
    print("\nStep 4: Saving JSON files...")
    for name in CELL_LINES:
        if name in expr_results:
            path = os.path.join(OUTPUT_DIR, f'depmap_expression_{name}.json')
            with open(path, 'w') as f:
                json.dump(expr_results[name], f, separators=(',', ':'))
            size_kb = os.path.getsize(path) / 1024
            print(f"  {path} ({size_kb:.0f} KB, {len(expr_results[name])} genes)")

        if name in crispr_results:
            path = os.path.join(OUTPUT_DIR, f'depmap_crispr_{name}.json')
            with open(path, 'w') as f:
                json.dump(crispr_results[name], f, separators=(',', ':'))
            size_kb = os.path.getsize(path) / 1024
            print(f"  {path} ({size_kb:.0f} KB, {len(crispr_results[name])} genes)")

    print("\nDone! Created 10 JSON files in web_data/")


if __name__ == '__main__':
    main()
