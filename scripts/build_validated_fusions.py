#!/usr/bin/env python3
"""
Build web_data/validated_fusions.json — the curated/validated fusion set used by
the mutational-analysis "Fusion" axis.

Source: web_data/clinical_fusions.json (curated driver fusions, validated per cell
line by tissue match + partner expression z-score + partner CRISPR dependency
z-score; each call carries a tier in {high, medium, low}). This is the same data
the Cell Line Browser shows under "Clinically relevant fusions".

We keep only HIGH + MEDIUM tier calls (drop the handful of low-tier ones) and key
by the SPECIFIC FUSION (e.g. "BCR-ABL1") — the analysis lets the user pick a named
driver fusion, not a gene. Presence is binary per fusion (a line carries it or it
does not). Output MIRRORS translocations.json so the fusion-axis code consumes it
unchanged, except the keys are fusion names instead of gene symbols:

{
  "genes": [...],                       # fusion names (e.g. "BCR-ABL1")
  "geneCounts": {fusion: N_lines},
  "geneData": {
    "BCR-ABL1": {
      "translocations": {cellLine: 1, ...},     # carriers (binary)
      "partners":       {cellLine: [geneA, geneB]},  # participating genes (hover)
      "counts": {"0": N_wt, "1": N, "2": 0},
      "total_translocated": N
    }
  }
}

Unlike process_translocations.py there is NO minimum-cell-line filter — every
validated fusion is curated and kept, however rare.

Usage:
    python build_validated_fusions.py
"""

import json
import os
from collections import defaultdict

ALLOWED_TIERS = {"high", "medium"}


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    metadata_file = os.path.join(project_dir, "web_data", "metadata.json")
    clinical_file = os.path.join(project_dir, "web_data", "clinical_fusions.json")
    output_file = os.path.join(project_dir, "web_data", "validated_fusions.json")

    with open(metadata_file) as f:
        cohort = set(json.load(f)["cellLines"])
    print(f"Loaded {len(cohort)} cell lines from metadata.json")

    with open(clinical_file) as f:
        clinical = json.load(f)
    fusion_data = clinical["fusionData"]
    print(f"{len(fusion_data)} curated fusions in clinical_fusions.json")

    # Key by the SPECIFIC FUSION (e.g. "BCR-ABL1"), not by participating gene —
    # the analysis axis lets the user pick a named driver fusion and splits cells
    # into carriers (1) vs non-carriers (0). Presence is binary per fusion.
    fusion_lines = defaultdict(set)  # fusion name -> set of cell lines (high/med tier)
    fusion_genes = {}                # fusion name -> participating gene list (for hover)
    n_calls = 0
    for fusion_name, fd in fusion_data.items():
        fusion_genes[fusion_name] = [g.strip() for g in fusion_name.split("-") if g.strip()]
        for cell_line, info in fd.get("cellLines", {}).items():
            if info.get("tier") not in ALLOWED_TIERS:
                continue
            if cell_line not in cohort:
                continue
            n_calls += 1
            fusion_lines[fusion_name].add(cell_line)

    print(f"{n_calls} high/medium-tier per-line calls kept")

    names_sorted = sorted(fusion_lines.keys())
    gene_counts = {}
    gene_data = {}
    for name in names_sorted:
        cls = fusion_lines[name]
        total = len(cls)
        gene_counts[name] = total
        # Binary presence: a line carries this fusion (1) or it does not (0).
        translocations = {cl: 1 for cl in sorted(cls)}
        partners = {cl: fusion_genes.get(name, []) for cl in sorted(cls)}
        gene_data[name] = {
            "translocations": translocations,
            "partners": partners,
            "counts": {
                "0": len(cohort) - total,
                "1": total,
                "2": 0,
            },
            "total_translocated": total,
        }

    output = {"genes": names_sorted, "geneCounts": gene_counts, "geneData": gene_data}
    with open(output_file, "w") as f:
        json.dump(output, f, separators=(",", ":"))

    print(f"Written {output_file}")
    print(f"  {len(names_sorted)} fusions")
    print(f"  Top genes by line count: {sorted(gene_counts.items(), key=lambda x: -x[1])[:12]}")


if __name__ == "__main__":
    main()
