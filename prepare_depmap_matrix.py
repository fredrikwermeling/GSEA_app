#!/usr/bin/env python3
"""
Build the DepMap matrices that Enrich reads one cell line at a time.

Output (web_data/):
  depmap_crispr.bin      int16, row-major, one row per cell line: Chronos gene effect x 1000
  depmap_expression.bin  int16, row-major, one row per cell line: log2FC vs the
                         median of all cell lines x 1000 (matrix is log2(TPM+1))
  depmap_index.json      gene lists, cell line list (id, name, lineage, disease),
                         row sizes and scale, so the app can fetch a single row
                         with an HTTP Range request (~37 KB) instead of a file per line.

Inputs in web_data/.cache (copy or symlink the DepMap release files there):
  CRISPRGeneEffect.csv, OmicsExpressionProteinCodingGenesTPMLogp1.csv, Model26Q1.csv
"""
import json, os, sys
import numpy as np
import pandas as pd

CACHE = 'web_data/.cache'
OUT = 'web_data'
RELEASE = '26Q1'
# Cell lines left out on purpose. HeLa: the lab does not distribute HeLa-derived
# data (Henrietta Lacks family agreement), whatever the DepMap release contains.
EXCLUDE = {'ACH-001086'}
SCALE = 1000
NA = -32768

def clean_gene(col):
    return col.split('(')[0].strip() if '(' in col else col.strip()

def load_models():
    m = pd.read_csv(os.path.join(CACHE, 'Model26Q1.csv'), dtype=str).fillna('')
    return {r['ModelID']: {'name': r['CellLineName'] or r['StrippedCellLineName'] or r['ModelID'],
                           'stripped': r['StrippedCellLineName'],
                           'lineage': r['OncotreeLineage'], 'disease': r['OncotreePrimaryDisease']}
            for _, r in m.iterrows()}

def write_matrix(df, models, name, transform):
    df = df.select_dtypes(include=[np.number])
    genes = [clean_gene(c) for c in df.columns]
    # keep the first column for a duplicated symbol
    keep = ~pd.Index(genes).duplicated()
    df = df.loc[:, keep]; genes = [g for g, k in zip(genes, keep) if k]
    df = df[~df.index.duplicated(keep='first')]
    ids = [i for i in df.index if i in models and i not in EXCLUDE]
    df = df.loc[ids]
    df = df.iloc[np.argsort([models[i]['name'].upper() for i in ids])]
    values = transform(df)
    arr = np.where(np.isnan(values), NA, np.clip(np.round(values * SCALE), -32767, 32767)).astype('<i2')
    path = os.path.join(OUT, f'depmap_{name}.bin')
    arr.tofile(path)
    print(f'  {path}: {arr.shape[0]} cell lines x {arr.shape[1]} genes, {os.path.getsize(path)/1e6:.1f} MB')
    return {'genes': genes,
            'cellLines': [{'id': i, 'name': models[i]['name'], 'lineage': models[i]['lineage'], 'disease': models[i]['disease']} for i in df.index],
            'nGenes': len(genes), 'nCellLines': len(df.index), 'rowBytes': len(genes) * 2}

def main():
    models = load_models()
    print('CRISPR...')
    crispr = pd.read_csv(os.path.join(CACHE, 'CRISPRGeneEffect.csv'), index_col=0)
    idx_c = write_matrix(crispr, models, 'crispr', lambda d: d.values.astype(float))
    del crispr
    print('Expression...')
    raw = pd.read_csv(os.path.join(CACHE, 'OmicsExpressionProteinCodingGenesTPMLogp1.csv'), index_col=0)
    # The expression file has one row per sequencing profile, so a model can
    # appear several times. Keep DepMap's default profile per model (or the
    # first row when the flag column is missing) so each cell line counts once
    # in the median and appears once in the list.
    if 'IsDefaultEntryForModel' in raw.columns:
        raw = raw[raw['IsDefaultEntryForModel'].astype(str).str.lower().isin(['yes', 'true', '1'])]
    if 'ModelID' in raw.columns:
        raw = raw.set_index('ModelID')
    raw = raw[~raw.index.duplicated(keep='first')]
    expr = raw.select_dtypes(include=[np.number])
    # median across all cell lines, so a row is log2FC vs the panel
    idx_e = write_matrix(expr, models, 'expression', lambda d: (d - d.median(axis=0)).values.astype(float))
    index = {'release': RELEASE, 'scale': SCALE, 'na': NA,
             'crispr': dict(idx_c, metric='Chronos_score', label='CRISPR screen, Chronos gene effect (negative = the cell line depends on the gene)'),
             'expression': dict(idx_e, metric='log2FC_vs_median', label='expression as log2 fold change vs the median of all DepMap cell lines')}
    with open(os.path.join(OUT, 'depmap_index.json'), 'w') as f:
        json.dump(index, f, separators=(',', ':'))
    print(f'  web_data/depmap_index.json: {os.path.getsize(os.path.join(OUT, "depmap_index.json"))/1e6:.1f} MB')

if __name__ == '__main__':
    main()
