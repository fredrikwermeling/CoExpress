#!/usr/bin/env python3
"""
Generate web_data/inferred_subtypes.json from DepMap's
OmicsInferredMolecularSubtypes.csv.

This file is DepMap's own per-cell-line annotation of clinically relevant
features inferred from integrated omics:
  * Specific hotspot variants (KRAS p.G12D, BRAF p.V600E, etc.)
  * Named driver fusions called by their pipeline
  * LoF (loss of function) calls for tumor suppressors — combines CN, mutation
    and expression, so it catches deletion-driven losses (e.g. CDKN2A homo-del)
    that a damaging-mutation matrix alone would miss
  * MSI status

Cells values: 'True', 'False', or '' (NaN — no omics data for that cell line).

Output (compact JSON):
{
  "schemaVersion": "1.0",
  "hotspots": [...column names like "KRAS p.G12D"...],
  "fusions":  [...column names like "BCR-ABL1"...],
  "lofGenes": [...gene names like "CDKN2A"...],
  "byCellLine": {
    "ACH-000322": {
      "hotspots": ["BRAF p.V600E"],
      "fusions":  [],
      "lof":      ["CDKN2A", "MTAP"],
      "msi":      false
    }, ...
  }
}

Usage:
  python scripts/process_inferred_subtypes.py <path_to_OmicsInferredMolecularSubtypes.csv>
"""

import csv
import json
import os
import sys
from collections import defaultdict


# Column buckets — keep in sync with the file header. The ordering here
# determines display order in the cell line browser.
HOTSPOT_COLS = [
    "KRAS p.G12D", "KRAS p.G12C", "BRAF p.V600E", "EGFR p.L858R", "JAK2 p.V617F",
    "KRAS p.G12", "KRAS p.G13", "KRAS p.Q61",
    "NRAS p.G12", "NRAS p.G13", "NRAS p.Q61",
    "HRAS p.G12", "HRAS p.G13", "HRAS p.Q61",
    "PIK3CA p.E542", "PIK3CA p.E545", "PIK3CA p.H1047",
    "ALK Hotspot", "EGFR exon 19 del",
]

FUSION_COLS = [
    "EWSR1-FLI1", "EWSR1-ERG", "EWSR1-FEV", "CIC-DUX4", "LMO2-STAG2",
    "ETV6-RUNX1", "TCF3-PBX1", "BCR-ABL1", "KMT2A Fusions", "PAX-FOXO1",
    "MEF2D Fusions",
]

LOF_COLS = [
    "RB1_LoF", "TP53_LoF", "PTEN_LoF", "NF1_LoF",
    "CDKN2A_LoF", "VHL_LoF", "MTAP_LoF", "APC_LoF",
]

MSI_COL = "MSI"


def parse_bool(s):
    """Tri-state: 'True' -> True, 'False' -> False, '' -> None (no omics)."""
    if s is None:
        return None
    t = s.strip()
    if t == "":
        return None
    return t.lower() == "true"


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: python process_inferred_subtypes.py <OmicsInferredMolecularSubtypes.csv>")
    csv_path = sys.argv[1]

    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    web = os.path.join(project_dir, "web_data")

    # Load valid cell lines (only emit data for lines we actually have GE for)
    with open(os.path.join(web, "metadata.json")) as f:
        ge_meta = json.load(f)
    valid_cell_lines = set(ge_meta["cellLines"])

    by_cl = {}
    n_total = n_kept = 0
    column_counts = defaultdict(int)

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        cols = reader.fieldnames or []
        # Sanity-check columns
        for c in HOTSPOT_COLS + FUSION_COLS + LOF_COLS + [MSI_COL]:
            if c not in cols:
                print(f"  warning: column not found in CSV: {c!r}")

        for row in reader:
            n_total += 1
            cl = (row.get("ModelID") or "").strip()
            if not cl or cl not in valid_cell_lines:
                continue
            n_kept += 1

            entry = {}

            hs = [c for c in HOTSPOT_COLS if parse_bool(row.get(c)) is True]
            if hs:
                entry["hotspots"] = hs
                for c in hs:
                    column_counts[c] += 1

            fus = [c for c in FUSION_COLS if parse_bool(row.get(c)) is True]
            if fus:
                entry["fusions"] = fus
                for c in fus:
                    column_counts[c] += 1

            lof = [c.replace("_LoF", "") for c in LOF_COLS if parse_bool(row.get(c)) is True]
            if lof:
                entry["lof"] = lof
                for c in lof:
                    column_counts[f"{c}_LoF"] += 1

            msi_val = parse_bool(row.get(MSI_COL))
            if msi_val is True:
                entry["msi"] = True
                column_counts["MSI"] += 1
            elif msi_val is False:
                # Keep explicit False so the UI can distinguish "not MSI" from
                # "no omics data". Saves bytes by omitting None.
                entry["msi"] = False

            if entry:  # don't store empty entries
                by_cl[cl] = entry

    output = {
        "_description": (
            "DepMap inferred molecular subtypes (per-cell-line specific hotspots, "
            "named driver fusions, integrated LoF for tumor suppressors, MSI). "
            "Source: OmicsInferredMolecularSubtypes.csv. LoF combines CN < 0.3, "
            "Likely-LoF mutation with AF > 0.5, or expression < 0.1 log-TPM — so "
            "it catches deletion-driven losses (e.g. CDKN2A homozygous deletion) "
            "that a damaging-mutation matrix alone would miss."
        ),
        "schemaVersion": "1.0",
        "hotspots": HOTSPOT_COLS,
        "fusions": FUSION_COLS,
        "lofGenes": [c.replace("_LoF", "") for c in LOF_COLS],
        "byCellLine": by_cl,
    }

    out_path = os.path.join(web, "inferred_subtypes.json")
    with open(out_path, "w") as f:
        json.dump(output, f, separators=(",", ":"))

    print(f"\nWritten {out_path}")
    print(f"  {n_total} rows seen, {n_kept} kept (cell line in metadata)")
    print(f"  {len(by_cl)} cell lines with at least one feature")

    # Cross-check vs clinical_fusions.json — report agreement
    cf_path = os.path.join(web, "clinical_fusions.json")
    if os.path.exists(cf_path):
        with open(cf_path) as f:
            cf = json.load(f)
        print(f"\nCross-check vs clinical_fusions.json:")
        # Map DepMap fusion column → our canonical fusion(s)
        overlap_map = {
            "BCR-ABL1":     ["BCR-ABL1"],
            "EWSR1-FLI1":   ["EWSR1-FLI1"],
            "EWSR1-ERG":    ["EWSR1-ERG"],
            "ETV6-RUNX1":   ["ETV6-RUNX1"],
            "TCF3-PBX1":    ["TCF3-PBX1"],
            "CIC-DUX4":     ["CIC-DUX4"],
            "PAX-FOXO1":    ["PAX3-FOXO1", "PAX7-FOXO1"],
            "KMT2A Fusions": [k for k in cf.get("fusionData", {}).keys() if k.startswith("KMT2A-")],
        }
        for depmap_name, our_names in overlap_map.items():
            depmap_set = {cl for cl, e in by_cl.items() if depmap_name in e.get("fusions", [])}
            our_set = set()
            for n in our_names:
                our_set.update((cf.get("fusionData", {}).get(n, {}).get("cellLines") or {}).keys())
            inter = depmap_set & our_set
            only_dep = depmap_set - our_set
            only_us = our_set - depmap_set
            print(f"  {depmap_name:20s}  DepMap={len(depmap_set):3d}  ours={len(our_set):3d}  agree={len(inter):3d}  DepMap-only={len(only_dep)}  ours-only={len(only_us)}")

    print(f"\nTop column counts (cell lines positive):")
    for c, n in sorted(column_counts.items(), key=lambda x: -x[1])[:25]:
        print(f"  {c:30s}  n={n}")


if __name__ == "__main__":
    main()
