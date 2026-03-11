#!/usr/bin/env python3
"""Convert MSigDB GMT file to compact JSON format for the GSEA web app.

Usage:
    python3 convert_gmt.py input.gmt output.json

The GMT format is tab-separated:
    GENE_SET_NAME<tab>DESCRIPTION<tab>GENE1<tab>GENE2<tab>...

The output JSON format is:
    { "GENE_SET_NAME": ["GENE1", "GENE2", ...], ... }
"""

import json
import sys
import os


def convert_gmt_to_json(gmt_path, json_path):
    gene_sets = {}
    with open(gmt_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                print(f"  Warning: line {line_num} has fewer than 3 columns, skipping")
                continue
            name = parts[0].strip()
            # parts[1] is the description (usually a URL), skip it
            genes = [g.strip().upper() for g in parts[2:] if g.strip()]
            if genes:
                gene_sets[name] = genes

    with open(json_path, 'w') as f:
        json.dump(gene_sets, f)

    size_mb = os.path.getsize(json_path) / (1024 * 1024)
    print(f"Converted {len(gene_sets)} gene sets -> {json_path} ({size_mb:.1f} MB)")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python3 convert_gmt.py <input.gmt> <output.json>")
        sys.exit(1)
    convert_gmt_to_json(sys.argv[1], sys.argv[2])
