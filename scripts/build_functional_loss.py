#!/usr/bin/env python3
"""
Build web_data/functional_loss.json — the validated tumor-suppressor functional-loss
set used by the mutational-analysis "Functional loss" axis (replacing the noisy
predicted-damaging matrix as a primary axis).

Source: web_data/inferred_subtypes.json, field byCellLine[cl].lof — DepMap's
integrated functional-loss call per tumor suppressor, where a gene is "lost" if any
of: deep CN loss (<= 0.3x diploid), a loss-of-function mutation with AF > 0.5, or
near-zero expression (< 0.1 log-TPM). This catches deletion-driven losses the
damaging-mutation matrix misses. It is the same data the Cell Line Browser shows
under "Functional loss". Covers 8 TSGs: APC, CDKN2A, MTAP, NF1, PTEN, RB1, TP53, VHL.

Every cell line in the CRISPR cohort has an inferred-subtype entry, so absence of a
gene from a line's lof list is a true wild-type (functionally intact) call — no
coverage mask needed.

Output MIRRORS damaging_mutations.json exactly so the existing damaging-axis code
consumes it unchanged:

{
  "genes": [...],                  # sorted by line count descending
  "geneCounts": {gene: N},
  "geneData": {
    GENE: {
      "mutations": {cellLine: 1, ...},   # lines with functional loss of GENE
      "counts": {"0": N_wt, "1": N},
      "total_mutated": N
    }
  }
}

Usage:
    python build_functional_loss.py
"""

import json
import os
from collections import defaultdict

# Genes that are deletion/co-deletion markers rather than classic point-mutated TSGs
# are still included — the source already curates the set to 8 meaningful genes.


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    metadata_file = os.path.join(project_dir, "web_data", "metadata.json")
    inferred_file = os.path.join(project_dir, "web_data", "inferred_subtypes.json")
    output_file = os.path.join(project_dir, "web_data", "functional_loss.json")

    with open(metadata_file) as f:
        cohort = list(json.load(f)["cellLines"])
    cohort_set = set(cohort)
    print(f"Loaded {len(cohort)} cell lines from metadata.json")

    with open(inferred_file) as f:
        inferred = json.load(f)
    by_cell_line = inferred["byCellLine"]
    n_with_inferred = sum(1 for cl in cohort if cl in by_cell_line)
    print(f"{n_with_inferred}/{len(cohort)} cohort lines have inferred-subtype data")

    gene_mutations = defaultdict(dict)  # gene -> {cell_line: 1}
    for cl in cohort:
        lof = (by_cell_line.get(cl) or {}).get("lof") or []
        for gene in lof:
            gene_mutations[gene][cl] = 1

    gene_counts = {}
    gene_data = {}
    n_cohort = len(cohort_set)
    for gene, mutations in gene_mutations.items():
        n = len(mutations)
        gene_counts[gene] = n
        gene_data[gene] = {
            "mutations": mutations,
            "counts": {"0": n_cohort - n, "1": n},
            "total_mutated": n,
        }

    genes_sorted = sorted(gene_counts.keys(), key=lambda g: -gene_counts[g])
    output = {"genes": genes_sorted, "geneCounts": gene_counts, "geneData": gene_data}
    with open(output_file, "w") as f:
        json.dump(output, f, separators=(",", ":"))

    print(f"Written {output_file}")
    print(f"  {len(genes_sorted)} genes")
    print(f"  Counts: {[(g, gene_counts[g]) for g in genes_sorted]}")


if __name__ == "__main__":
    main()
