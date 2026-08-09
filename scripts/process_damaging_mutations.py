#!/usr/bin/env python3
"""
Process DepMap OmicsSomaticMutationsMatrixDamaging.csv into damaging_mutations.json.

This is a binary matrix (0/1) of damaging somatic mutations per cell line per gene.
Unlike hotspot mutations (mutations.json with 49 genes), this covers ~12K genes.

Output format matches mutations.json structure:
{
  "genes": [...],           # sorted by mutation count descending
  "geneCounts": {gene: N},  # number of mutated cell lines
  "geneData": {
    gene: {
      "mutations": {cellLine: 1, ...},  # only mutated cell lines
      "counts": {"0": N_wt, "1": N_mut},
      "total_mutated": N
    }
  }
}

Usage:
    python process_damaging_mutations.py [path_to_csv]

    If no path given, uses scripts/OmicsSomaticMutationsMatrixDamaging.csv
"""

import csv
import json
import sys
import os
import re

MIN_CELL_LINES = 3  # Minimum mutated cell lines for a gene to be included


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)

    if len(sys.argv) >= 2:
        csv_file = sys.argv[1]
    else:
        csv_file = os.path.join(script_dir, "OmicsSomaticMutationsMatrixDamaging.csv")

    metadata_file = os.path.join(project_dir, "web_data", "metadata.json")
    output_file = os.path.join(project_dir, "web_data", "damaging_mutations.json")

    # Load valid cell lines from metadata
    with open(metadata_file) as f:
        metadata = json.load(f)
    valid_cell_lines = set(metadata["cellLines"])
    print(f"Loaded {len(valid_cell_lines)} cell lines from metadata.json")

    # Read CSV header to get gene names
    with open(csv_file) as f:
        reader = csv.reader(f)
        header = next(reader)

    # 26Q1 switched these matrices to one row per PROFILE with metadata
    # columns (SequencingID, ModelID, IsDefaultEntryForModel, ...); older
    # releases had ModelID as the unnamed first column. Handle both.
    model_col = header.index("ModelID") if "ModelID" in header else 0
    default_col = header.index("IsDefaultEntryForModel") if "IsDefaultEntryForModel" in header else None

    # Parse gene names from header (format: "GENE (EntrezID)"); metadata
    # columns don't match the pattern and are skipped.
    gene_cols = []  # (col_index, gene_name)
    for i, col in enumerate(header):
        match = re.match(r'^(.+?)\s*\(\d+\)$', col.strip())
        if match:
            gene_cols.append((i, match.group(1).strip()))

    print(f"CSV has {len(gene_cols)} genes")

    # Read data — collect mutations per gene
    # gene_name -> {cell_line: 1}
    gene_mutations = {name: {} for _, name in gene_cols}

    n_valid = 0
    n_total = 0
    with open(csv_file) as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            if default_col is not None and row[default_col].strip() != 'Yes':
                continue
            n_total += 1
            model_id = row[model_col].strip()
            if model_id not in valid_cell_lines:
                continue
            n_valid += 1
            for col_idx, gene_name in gene_cols:
                if col_idx < len(row):
                    v = row[col_idx].strip()
                    if v and v not in ('0', '0.0'):
                        try:
                            if float(v) > 0:
                                gene_mutations[gene_name][model_id] = 1
                        except ValueError:
                            pass

    print(f"Cell lines: {n_total} total, {n_valid} in metadata")

    # Filter by minimum cell lines and build output
    n_valid_cl = n_valid  # total valid cell lines for WT count
    gene_counts = {}
    gene_data = {}

    for gene_name, mutations in gene_mutations.items():
        n_mut = len(mutations)
        if n_mut < MIN_CELL_LINES:
            continue
        gene_counts[gene_name] = n_mut
        gene_data[gene_name] = {
            "mutations": mutations,
            "counts": {
                "0": n_valid_cl - n_mut,
                "1": n_mut,
            },
            "total_mutated": n_mut,
        }

    # Sort genes by mutation count descending
    genes_sorted = sorted(gene_counts.keys(), key=lambda g: -gene_counts[g])

    print(f"Genes with >= {MIN_CELL_LINES} mutated cell lines: {len(genes_sorted)}")
    print(f"Top 10: {[(g, gene_counts[g]) for g in genes_sorted[:10]]}")

    output = {
        "genes": genes_sorted,
        "geneCounts": gene_counts,
        "geneData": gene_data,
    }

    with open(output_file, "w") as f:
        json.dump(output, f, separators=(",", ":"))

    size_mb = os.path.getsize(output_file) / (1024 * 1024)
    print(f"Written {output_file} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
