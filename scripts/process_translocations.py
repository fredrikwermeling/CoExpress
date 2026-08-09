#!/usr/bin/env python3
"""
Process DepMap OmicsFusionFiltered.csv into translocations.json for Correlate V2.

Reads Arriba-pipeline fusion data, filters to cell lines present in metadata.json,
indexes by both gene1 and gene2, counts distinct fusion partners per cell line per gene,
and outputs a JSON file matching the mutations.json structure.

Usage:
    python process_translocations.py <path_to_OmicsFusionFiltered.csv>

Output: web_data/translocations.json
"""

import csv
import json
import sys
import os
from collections import defaultdict

MIN_CELL_LINES = 3  # Minimum cell lines with fusions for a gene to be included


def main():
    if len(sys.argv) < 2:
        print("Usage: python process_translocations.py <OmicsFusionFiltered.csv>")
        sys.exit(1)

    fusion_file = sys.argv[1]
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    metadata_file = os.path.join(project_dir, "web_data", "metadata.json")
    output_file = os.path.join(project_dir, "web_data", "translocations.json")

    # Load valid cell lines from metadata
    with open(metadata_file) as f:
        metadata = json.load(f)
    valid_cell_lines = set(metadata["cellLines"])
    print(f"Loaded {len(valid_cell_lines)} cell lines from metadata.json")

    # Read fusion data
    # Expected columns: ModelID (ACH-...), FusionName (e.g. BCR--ABL1),
    # LeftGene, RightGene (or Gene1, Gene2)
    # Try to detect column names from header
    fusions = []
    with open(fusion_file, newline="") as f:
        reader = csv.DictReader(f)
        cols = reader.fieldnames
        print(f"Columns found: {cols}")

        # Detect column names
        model_col = None
        for c in cols:
            if c in ("ModelID", "DepMap_ID", "CellLineID", "ModelId"):
                model_col = c
                break
        if model_col is None:
            # Try first column
            model_col = cols[0]
            print(f"Warning: using first column '{model_col}' as model ID")

        # Detect gene columns
        gene1_col = gene2_col = fusion_name_col = None
        for c in cols:
            cl = c.lower().replace("_", "").replace(" ", "")
            if cl in ("leftgene", "gene1", "genea", "leftgenename"):
                gene1_col = c
            elif cl in ("rightgene", "gene2", "geneb", "rightgenename"):
                gene2_col = c
            elif cl in ("fusionname", "fusion", "name"):
                fusion_name_col = c

        # If gene columns not found, try to parse from FusionName
        if gene1_col is None or gene2_col is None:
            if fusion_name_col:
                print(f"Gene columns not found directly, will parse from '{fusion_name_col}'")
            else:
                print("Error: Cannot find gene columns or fusion name column")
                sys.exit(1)

        for row in reader:
            model_id = row[model_col].strip()
            if model_id not in valid_cell_lines:
                continue

            if gene1_col and gene2_col:
                g1 = row[gene1_col].strip()
                g2 = row[gene2_col].strip()
                # Strip Ensembl IDs in parentheses, e.g. "BCR (ENSG00000186716.20)" -> "BCR"
                if " (" in g1:
                    g1 = g1.split(" (")[0].strip()
                if " (" in g2:
                    g2 = g2.split(" (")[0].strip()
            elif fusion_name_col:
                fname = row[fusion_name_col].strip()
                # Parse "BCR--ABL1" or "BCR-ABL1" format
                if "--" in fname:
                    parts = fname.split("--", 1)
                elif "-" in fname:
                    parts = fname.split("-", 1)
                else:
                    continue
                g1 = parts[0].strip()
                g2 = parts[1].strip() if len(parts) > 1 else ""
            else:
                continue

            if not g1 or not g2:
                continue

            fusions.append((model_id, g1, g2))

    print(f"Found {len(fusions)} fusion events in valid cell lines")

    # Build gene -> cell line -> set of partners
    # Index under both genes
    gene_partners = defaultdict(lambda: defaultdict(set))

    for model_id, g1, g2 in fusions:
        gene_partners[g1][model_id].add(g2)
        gene_partners[g2][model_id].add(g1)

    # Filter genes by minimum cell line threshold
    filtered_genes = {}
    for gene, cell_data in gene_partners.items():
        n_cell_lines = len(cell_data)
        if n_cell_lines >= MIN_CELL_LINES:
            filtered_genes[gene] = cell_data

    print(f"Genes with fusions in >= {MIN_CELL_LINES} cell lines: {len(filtered_genes)}")

    # Build output structure
    genes_sorted = sorted(filtered_genes.keys())
    gene_counts = {}
    gene_data = {}

    for gene in genes_sorted:
        cell_data = filtered_genes[gene]
        total_translocated = len(cell_data)
        gene_counts[gene] = total_translocated

        # Build translocations map (cell line -> number of distinct partners)
        translocations = {}
        partners = {}
        counts = defaultdict(int)

        for cell_line, partner_set in cell_data.items():
            n_partners = len(partner_set)
            translocations[cell_line] = n_partners
            partners[cell_line] = sorted(partner_set)
            bucket = str(min(n_partners, 2))  # 1 or 2+
            counts[bucket] += 1

        # Count WT (cell lines without fusions)
        n_wt = len(valid_cell_lines) - total_translocated
        counts["0"] = n_wt

        gene_data[gene] = {
            "translocations": translocations,
            "partners": partners,
            "counts": {
                "0": counts["0"],
                "1": counts.get("1", 0),
                "2": counts.get("2", 0),
            },
            "total_translocated": total_translocated,
        }

    output = {
        "genes": genes_sorted,
        "geneCounts": gene_counts,
        "geneData": gene_data,
    }

    with open(output_file, "w") as f:
        json.dump(output, f, separators=(",", ":"))

    print(f"Written {output_file}")
    print(f"  {len(genes_sorted)} genes")
    print(f"  Top genes by fusion count: {sorted(gene_counts.items(), key=lambda x: -x[1])[:10]}")


if __name__ == "__main__":
    main()
