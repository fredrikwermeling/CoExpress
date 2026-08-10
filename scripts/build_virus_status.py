#!/usr/bin/env python3
"""
Build web_data/virus_status.json from Cellosaurus.

Cellosaurus records what a line was transformed by in a structured comment:

    CC   Transformant: NCBI_TaxID; 10376; Epstein-Barr virus (EBV).

That is a CONFIRMED, curated statement about a line. It is not a survey of
infection: a line with no Transformant comment has simply not been recorded as
virus-transformed, which is not the same as being virus-free. Every label the
app shows is worded accordingly.

Chemical and physical transformants (MNNG, X-Ray, ...) are ignored; only viral
agents are kept, matched by their NCBI taxonomy id.

Usage:
    python3 scripts/build_virus_status.py /path/to/cellosaurus.txt

Output: web_data/virus_status.json
    {
      "_description": ..., "source": ..., "retrieved": "YYYY-MM-DD",
      "agents": { "EBV": {"label": ..., "taxid": ...}, ... },
      "byCellLine": { "ACH-000001": [ {"agent": "EBV", "name": ..., "note": ...} ] },
      "counts": { "EBV": 42, ... }
    }
"""

import json
import os
import re
import sys
from datetime import date

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
WEB = os.path.join(ROOT, "web_data")

# Viral transformants worth surfacing, by NCBI taxonomy id. Chemical and
# physical agents in the same field (MNNG, MCA, X-Ray) are deliberately absent.
AGENTS = {
    "10376": ("EBV", "Epstein-Barr virus"),
    "333761": ("HPV18", "Human papillomavirus type 18"),
    "333760": ("HPV16", "Human papillomavirus type 16"),
    "10566": ("HPV", "Human papillomavirus"),
    "1891762": ("SV40", "Simian virus 40"),
    "10633": ("SV40", "Simian virus 40"),
    "11908": ("HTLV-1", "Human T-lymphotropic virus 1"),
    "37296": ("KSHV", "Kaposi sarcoma-associated herpesvirus (HHV-8)"),
    "10407": ("HBV", "Hepatitis B virus"),
    "11103": ("HCV", "Hepatitis C virus"),
    "11676": ("HIV-1", "Human immunodeficiency virus 1"),
    "28285": ("AdV5", "Human adenovirus C serotype 5"),
    "11801": ("MoMuLV", "Moloney murine leukemia virus"),
    "11788": ("AbMuLV", "Abelson murine leukemia virus"),
}

TRANSFORMANT = re.compile(r"^CC   Transformant: NCBI_TaxID; (\d+); (.+?)\.?$")


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/cellosaurus.txt"
    meta = json.load(open(os.path.join(WEB, "cellLineMetadata.json")))
    cohort = meta["cellLines"]
    # Cellosaurus is keyed by RRID; our cohort carries one per line.
    rrid_to_model = {}
    for cl in cohort:
        r = (meta.get("rrid") or {}).get(cl)
        if r:
            rrid_to_model[r.strip()] = cl

    by_line, counts = {}, {}
    cur_rrid, cur_hits = None, []

    def flush():
        if not cur_rrid or not cur_hits:
            return
        model = rrid_to_model.get(cur_rrid)
        if not model:
            return
        seen, keep = set(), []
        for h in cur_hits:
            if h["agent"] in seen:
                continue
            seen.add(h["agent"])
            keep.append(h)
            counts[h["agent"]] = counts.get(h["agent"], 0) + 1
        by_line[model] = keep

    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("AC   "):
                flush()
                cur_rrid, cur_hits = line[5:].strip(), []
            elif line.startswith("CC   Transformant:"):
                m = TRANSFORMANT.match(line)
                if not m:
                    continue
                taxid, text = m.group(1), m.group(2)
                if taxid not in AGENTS:
                    continue
                code, name = AGENTS[taxid]
                note = ""
                nm = re.search(r"\(Note=([^)]*)\)", text)
                if nm:
                    note = nm.group(1)
                cur_hits.append({"agent": code, "name": name, "note": note})
            elif line.startswith("//"):
                flush()
                cur_rrid, cur_hits = None, []
    flush()

    out = {
        "_description": (
            "Curated virus-transformation status from Cellosaurus (the "
            "Transformant annotation). A cell line listed here is CONFIRMED to "
            "have been transformed by that agent. A line NOT listed has no such "
            "curated record, which is not evidence that it is free of the virus: "
            "most lines have never been surveyed for infection."
        ),
        "source": "Cellosaurus (https://www.cellosaurus.org), Transformant annotation",
        "retrieved": date.today().isoformat(),
        "agents": {c: {"label": n, "taxid": t} for t, (c, n) in AGENTS.items()},
        "counts": counts,
        "byCellLine": by_line,
    }
    dest = os.path.join(WEB, "virus_status.json")
    json.dump(out, open(dest, "w"))
    print(f"{len(by_line)} of {len(cohort)} cohort lines have a curated transformant")
    for code, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"  {code:8s} {n}")
    print("->", dest)


if __name__ == "__main__":
    main()
