#!/usr/bin/env python3
"""
Generate web_data/global_signatures.json from DepMap's OmicsGlobalSignatures.csv.

Per-cell-line genome-wide features computed once for the panel:
  * Ploidy   — average chromosome dosage from PureCN
  * WGD      — whole genome doubling, binary
  * CIN      — chromosomal instability score
  * LoH      — loss of heterozygosity fraction
  * MSIScore — MSIsensor2 score (>=20 == MSI-high; the inferred subtypes
               file already exposes the binary call, this is the underlying float)
  * Aneuploidy — Ben-David 2021 aneuploidy score (integer)

Output: web_data/global_signatures.json keyed by ModelID. Filters to default
entries (IsDefaultEntryForModel == Yes) and to cell lines we have GE data for.

Usage:
  python scripts/process_global_signatures.py <path_to_OmicsGlobalSignatures.csv>
"""

import csv
import json
import os
import sys


FIELDS = ["MSIScore", "LoHFraction", "WGD", "CIN", "Ploidy", "Aneuploidy"]


def parse_float(s):
    if s is None:
        return None
    t = s.strip()
    if t == "" or t.lower() in {"nan", "na", "none"}:
        return None
    try:
        return float(t)
    except ValueError:
        return None


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: python process_global_signatures.py <OmicsGlobalSignatures.csv>")
    csv_path = sys.argv[1]

    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    web = os.path.join(project_dir, "web_data")

    with open(os.path.join(web, "metadata.json")) as f:
        valid_cell_lines = set(json.load(f)["cellLines"])

    by_cl = {}
    n_total = n_kept = 0
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            n_total += 1
            cl = (row.get("ModelID") or "").strip()
            if not cl or cl not in valid_cell_lines:
                continue
            # Each cell line may have multiple sequencings; only use the
            # one DepMap marks as the default for the model.
            if (row.get("IsDefaultEntryForModel") or "").strip().lower() != "yes":
                continue
            n_kept += 1
            entry = {}
            for f_name in FIELDS:
                v = parse_float(row.get(f_name))
                if v is None:
                    continue
                if f_name == "WGD":
                    entry[f_name] = bool(round(v))
                elif f_name == "Aneuploidy":
                    entry[f_name] = int(round(v))
                else:
                    entry[f_name] = round(v, 3)
            if entry:
                by_cl[cl] = entry

    output = {
        "_description": (
            "Per-cell-line genome-wide signatures from DepMap "
            "OmicsGlobalSignatures.csv (PureCN + MSIsensor2 + Ben-David "
            "aneuploidy score). Default sequencing per model only."
        ),
        "schemaVersion": "1.0",
        "fields": FIELDS,
        "byCellLine": by_cl,
    }

    out_path = os.path.join(web, "global_signatures.json")
    with open(out_path, "w") as f:
        json.dump(output, f, separators=(",", ":"))

    # Summary
    n_wgd = sum(1 for e in by_cl.values() if e.get("WGD"))
    n_msi = sum(1 for e in by_cl.values() if (e.get("MSIScore") or 0) >= 20)
    print(f"Written {out_path}")
    print(f"  {n_total} rows seen, {n_kept} kept (default entries in metadata)")
    print(f"  {len(by_cl)} cell lines with at least one signature")
    print(f"  WGD-positive: {n_wgd}")
    print(f"  MSI-high (score >= 20): {n_msi}")


if __name__ == "__main__":
    main()
