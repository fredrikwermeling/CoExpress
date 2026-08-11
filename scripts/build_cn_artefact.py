#!/usr/bin/env python3
"""
Build web_data/cn_artefact.json: how strongly each gene's CRISPR score tracks
its own copy number.

The best-known false signal in a CRISPR screen is the copy-number effect. In an
amplified region every cut is lethal whatever the gene does, so a gene sitting
in an amplification looks essential wherever that amplification occurs. Chronos
corrects much of this, but not all of it, and the residue is exactly the kind of
"strong hit" that wastes a follow-up experiment.

DepMap does not ship a list of artefact genes in the files this app is built
from, so the signal is measured here instead: Pearson r between a gene's effect
and its own relative copy number, across every cell line with both. A strongly
NEGATIVE r is the signature (more copies -> more cutting -> lower score).

Only genes at r <= -0.20 are written, which keeps the file small; everything
else is by definition unremarkable and needs no entry.

    python3 scripts/build_cn_artefact.py

Output: web_data/cn_artefact.json
    { "_description": ..., "threshold": -0.2, "computed": "YYYY-MM-DD",
      "genes": { "GENE": {"r": -0.42, "n": 1150}, ... } }
"""

import gzip
import json
import os
from datetime import date

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WEB = os.path.join(os.path.dirname(HERE), "web_data")

KEEP_AT_OR_BELOW = -0.20
MIN_PAIRS = 50


def load_matrix(bin_name, meta_name):
    meta = json.load(open(os.path.join(WEB, meta_name)))
    with gzip.open(os.path.join(WEB, bin_name), "rb") as fh:
        raw = np.frombuffer(fh.read(), dtype="<i2")
    n_g, n_c = meta["nGenes"], meta["nCellLines"]
    arr = raw[: n_g * n_c].reshape(n_g, n_c).astype(np.float32)
    arr[arr == meta["naValue"]] = np.nan
    arr /= meta["scaleFactor"]
    return meta, arr


def main():
    ge_meta, ge = load_matrix("geneEffects.bin.gz", "metadata.json")
    cn_meta, cn = load_matrix("cn.bin.gz", "cn_metadata.json")
    print(f"gene effect {ge.shape}, copy number {cn.shape}")

    # Cell lines present in both, and the column each sits in.
    cn_col = {cl: i for i, cl in enumerate(cn_meta["cellLines"])}
    shared = [(i, cn_col[cl]) for i, cl in enumerate(ge_meta["cellLines"]) if cl in cn_col]
    ge_cols = np.array([a for a, _ in shared])
    cn_cols = np.array([b for _, b in shared])
    print(f"{len(shared)} cell lines have both")

    cn_gene_row = {g: i for i, g in enumerate(cn_meta["genes"])}
    out = {}
    for gi, gene in enumerate(ge_meta["genes"]):
        ci = cn_gene_row.get(gene)
        if ci is None:
            continue
        x = ge[gi, ge_cols]
        y = cn[ci, cn_cols]
        ok = ~np.isnan(x) & ~np.isnan(y)
        n = int(ok.sum())
        if n < MIN_PAIRS:
            continue
        xa, ya = x[ok], y[ok]
        if xa.std() == 0 or ya.std() == 0:
            continue
        r = float(np.corrcoef(xa, ya)[0, 1])
        if r <= KEEP_AT_OR_BELOW:
            out[gene] = {"r": round(r, 3), "n": n}

    payload = {
        "_description": (
            "How strongly each gene's CRISPR score tracks its own relative copy number, "
            "as Pearson r across cell lines with both measured. A strongly negative value "
            "is the signature of the copy-number cutting artefact: in an amplified region "
            "every cut is lethal whatever the gene does, so the gene reads as essential "
            "simply for sitting in an amplification. Computed from the data in this app; "
            "DepMap does not publish such a list in these files. Only genes at or below "
            "the threshold are listed, so absence means the gene shows no such pattern."
        ),
        "threshold": KEEP_AT_OR_BELOW,
        "minPairs": MIN_PAIRS,
        "computed": date.today().isoformat(),
        "genes": out,
    }
    dest = os.path.join(WEB, "cn_artefact.json")
    json.dump(payload, open(dest, "w"))
    print(f"{len(out)} genes flagged at r <= {KEEP_AT_OR_BELOW}")
    worst = sorted(out.items(), key=lambda kv: kv[1]["r"])[:12]
    for g, v in worst:
        print(f"  {g:12s} r={v['r']:.3f}  n={v['n']}")
    print("->", dest)


if __name__ == "__main__":
    main()
