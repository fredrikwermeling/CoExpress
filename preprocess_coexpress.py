#!/usr/bin/env python3
"""
Preprocess expression data for CoExpress web app.
Converts DepMap OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
into web-ready binary format matching Correlate's data structure.
"""

import csv
import json
import gzip
import struct
import re
import shutil
import os
import numpy as np

# Paths
EXPRESSION_CSV = "/Users/fredrikwermeling/Documents/coexpress/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
MODEL_CSV = "/Users/fredrikwermeling/Documents/correlate app feb 2026 (färdig)/correlation app green listed/Model_25Q3.csv"
MUTATION_CSV = "/Users/fredrikwermeling/Documents/correlate app feb 2026 (färdig)/correlation app green listed/OmicsSomaticMutationsMatrixHotspot_25Q3.csv"
CORRELATE_WEB_DATA = "/Users/fredrikwermeling/Documents/correlate app feb 2026 (färdig)/correlation app green listed/web_data"
OUTPUT_DIR = "/Users/fredrikwermeling/Documents/coexpress-app/web_data"

SCALE_FACTOR = 3000  # Expression values range 0-9.75, so 9.75*3000=29250 fits well in Int16
NA_VALUE = -32768

def parse_gene_name(col_header):
    """Extract gene name from 'GENENAME (EntrezID)' format."""
    match = re.match(r'^(.+?)\s*\(\d+\)$', col_header)
    if match:
        return match.group(1).strip()
    return col_header.strip()

def load_model_metadata():
    """Load cell line metadata from Model.csv."""
    print("Loading Model metadata...")
    metadata = {}
    with open(MODEL_CSV, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            model_id = row.get('ModelID', '')
            if model_id:
                metadata[model_id] = {
                    'cellLineName': row.get('CellLineName', ''),
                    'strippedCellLineName': row.get('StrippedCellLineName', ''),
                    'lineage': row.get('OncotreeLineage', ''),
                    'primaryDisease': row.get('OncotreePrimaryDisease', ''),
                    'subtype': row.get('OncotreeSubtype', ''),
                }
    print(f"  Loaded metadata for {len(metadata)} models")
    return metadata

def load_mutations():
    """Load hotspot mutation data from OmicsSomaticMutationsMatrixHotspot."""
    print("Loading mutation data...")

    # Read header to get gene columns
    with open(MUTATION_CSV, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)

    # Find gene columns (format: "GENE (EntrezID)")
    meta_cols = {'', 'ModelID', 'ProfileID', 'PatientID'}
    gene_cols = []
    gene_col_indices = []
    for i, col in enumerate(header):
        if col not in meta_cols and re.match(r'^.+\s*\(\d+\)$', col):
            gene_cols.append(col)
            gene_col_indices.append(i)

    # Find ModelID column
    model_id_col = header.index('ModelID') if 'ModelID' in header else None
    if model_id_col is None:
        print("  WARNING: No ModelID column found in mutation CSV")
        return {'genes': [], 'geneCounts': {}, 'geneData': {}}

    # Read mutation data
    gene_data = {}
    for col in gene_cols:
        gene_name = parse_gene_name(col)
        gene_data[gene_name] = {
            'column': col,
            'mutations': {},
            'counts': {'0': 0, '1': 0, '2': 0}
        }

    with open(MUTATION_CSV, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            model_id = row[model_id_col]
            for col_name, col_idx in zip(gene_cols, gene_col_indices):
                gene_name = parse_gene_name(col_name)
                val = row[col_idx].strip() if col_idx < len(row) else '0'
                try:
                    val_int = int(float(val)) if val else 0
                except ValueError:
                    val_int = 0
                val_int = min(val_int, 2)  # Cap at 2
                gene_data[gene_name]['counts'][str(val_int)] = gene_data[gene_name]['counts'].get(str(val_int), 0) + 1
                if val_int > 0:
                    gene_data[gene_name]['mutations'][model_id] = val_int

    # Filter to genes with at least 3 mutations (matching Correlate)
    filtered_genes = {}
    for gene_name, data in gene_data.items():
        total_mutated = len(data['mutations'])
        if total_mutated >= 3:
            data['total_mutated'] = total_mutated
            filtered_genes[gene_name] = data

    # Sort by total mutations
    sorted_genes = sorted(filtered_genes.keys(), key=lambda g: filtered_genes[g]['total_mutated'], reverse=True)

    # Take top 49 most mutated (matching Correlate)
    top_genes = sorted_genes[:49]

    result = {
        'genes': top_genes,
        'geneCounts': {g: filtered_genes[g]['total_mutated'] for g in top_genes},
        'geneData': {g: filtered_genes[g] for g in top_genes}
    }

    print(f"  Found {len(top_genes)} hotspot mutation genes")
    return result

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Step 1: Load model metadata
    model_meta = load_model_metadata()

    # Step 2: Read expression data
    print("Reading expression CSV...")
    with open(EXPRESSION_CSV, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)

    # Parse column structure
    # Columns: ,SequencingID,ModelID,IsDefaultEntryForModel,ModelConditionID,IsDefaultEntryForMC,GENE1 (ID),...
    meta_col_names = header[:6]  # first 6 columns are metadata
    gene_columns = header[6:]
    gene_names = [parse_gene_name(col) for col in gene_columns]
    n_genes = len(gene_names)
    print(f"  Found {n_genes} genes")

    # Find column indices
    model_id_col = header.index('ModelID')
    is_default_col = header.index('IsDefaultEntryForModel')

    # Read data, filtering to default entries
    print("Reading expression values (filtering to IsDefaultEntryForModel=Yes)...")
    cell_lines = []
    expression_data = []  # list of rows, each row is a list of floats

    with open(EXPRESSION_CSV, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row_num, row in enumerate(reader):
            if row_num % 200 == 0:
                print(f"  Processing row {row_num}...")

            is_default = row[is_default_col].strip()
            if is_default != 'Yes':
                continue

            model_id = row[model_id_col].strip()
            cell_lines.append(model_id)

            # Parse expression values
            values = []
            for i in range(6, len(row)):
                val = row[i].strip()
                try:
                    values.append(float(val))
                except (ValueError, IndexError):
                    values.append(float('nan'))

            # Pad if necessary
            while len(values) < n_genes:
                values.append(float('nan'))

            expression_data.append(values)

    n_cell_lines = len(cell_lines)
    print(f"  Loaded {n_cell_lines} cell lines x {n_genes} genes")

    # Step 3: Convert to numpy and quantize to Int16
    print("Converting to Int16 binary matrix...")
    data_array = np.array(expression_data, dtype=np.float32)

    # Quantize: value * SCALE_FACTOR, NaN -> NA_VALUE
    int16_data = np.full(data_array.shape, NA_VALUE, dtype=np.int16)
    valid_mask = ~np.isnan(data_array)
    int16_data[valid_mask] = np.clip(
        np.round(data_array[valid_mask] * SCALE_FACTOR).astype(np.int32),
        -32767, 32767
    ).astype(np.int16)

    # Transpose to gene-major order: [nGenes x nCellLines]
    # (matching Correlate's format where data is stored as gene rows)
    int16_transposed = int16_data.T  # Now [nGenes x nCellLines]

    # Flatten and write as gzipped binary
    flat_data = int16_transposed.flatten()
    raw_bytes = flat_data.tobytes()

    output_bin = os.path.join(OUTPUT_DIR, 'expression.bin.gz')
    print(f"  Writing {output_bin}...")
    with gzip.open(output_bin, 'wb') as f:
        f.write(raw_bytes)

    uncompressed_mb = len(raw_bytes) / (1024 * 1024)
    compressed_mb = os.path.getsize(output_bin) / (1024 * 1024)
    print(f"  Uncompressed: {uncompressed_mb:.1f} MB, Compressed: {compressed_mb:.1f} MB")

    # Step 4: Create metadata.json
    print("Creating metadata.json...")

    # genesFull includes the original column names with EntrezID
    genes_full = gene_columns

    metadata = {
        'nGenes': n_genes,
        'nCellLines': n_cell_lines,
        'scaleFactor': SCALE_FACTOR,
        'naValue': NA_VALUE,
        'genes': gene_names,
        'genesFull': genes_full,
        'cellLines': cell_lines
    }

    with open(os.path.join(OUTPUT_DIR, 'metadata.json'), 'w') as f:
        json.dump(metadata, f)
    print(f"  {n_genes} genes, {n_cell_lines} cell lines")

    # Step 5: Create cellLineMetadata.json
    print("Creating cellLineMetadata.json...")

    cell_line_meta = {
        'cellLines': cell_lines,
        'cellLineName': {},
        'strippedCellLineName': {},
        'lineage': {},
        'primaryDisease': {},
        'subtype': {}
    }

    matched = 0
    for cl in cell_lines:
        if cl in model_meta:
            m = model_meta[cl]
            cell_line_meta['cellLineName'][cl] = m['cellLineName']
            cell_line_meta['strippedCellLineName'][cl] = m['strippedCellLineName']
            cell_line_meta['lineage'][cl] = m['lineage']
            cell_line_meta['primaryDisease'][cl] = m['primaryDisease']
            cell_line_meta['subtype'][cl] = m['subtype']
            matched += 1
        else:
            cell_line_meta['cellLineName'][cl] = cl
            cell_line_meta['strippedCellLineName'][cl] = cl
            cell_line_meta['lineage'][cl] = 'Unknown'
            cell_line_meta['primaryDisease'][cl] = 'Unknown'
            cell_line_meta['subtype'][cl] = ''

    print(f"  Matched {matched}/{n_cell_lines} cell lines to metadata")

    with open(os.path.join(OUTPUT_DIR, 'cellLineMetadata.json'), 'w') as f:
        json.dump(cell_line_meta, f)

    # Step 6: Create mutations.json
    print("Creating mutations.json...")
    mutations = load_mutations()
    with open(os.path.join(OUTPUT_DIR, 'mutations.json'), 'w') as f:
        json.dump(mutations, f)

    # Step 7: Copy orthologs.json and synonyms.json from Correlate
    print("Copying orthologs.json and synonyms.json from Correlate...")
    for fname in ['orthologs.json', 'synonyms.json']:
        src = os.path.join(CORRELATE_WEB_DATA, fname)
        dst = os.path.join(OUTPUT_DIR, fname)
        if os.path.exists(src):
            shutil.copy2(src, dst)
            print(f"  Copied {fname}")
        else:
            print(f"  WARNING: {fname} not found at {src}")

    # Step 8: Copy logo
    print("Copying logo...")
    logo_src = "/Users/fredrikwermeling/Documents/coexpress/coexpress logo 1.png"
    logo_dst = os.path.join(OUTPUT_DIR, "coexpress_logo.png")
    if os.path.exists(logo_src):
        shutil.copy2(logo_src, logo_dst)
        print(f"  Copied logo to {logo_dst}")
    else:
        print(f"  WARNING: Logo not found at {logo_src}")

    print("\n=== DONE ===")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"  expression.bin.gz: {compressed_mb:.1f} MB")
    print(f"  metadata.json: {n_genes} genes x {n_cell_lines} cell lines")
    print(f"  cellLineMetadata.json: {matched} cell lines matched")
    print(f"  mutations.json: {len(mutations['genes'])} hotspot genes")

if __name__ == '__main__':
    main()
