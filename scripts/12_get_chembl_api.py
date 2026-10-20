#!/usr/bin/env python
"""
NP-BioSynAct -- pull ChEMBL bioactivity for a pain / chronic-disease target panel
via the REST API. No 5.5GB download; minutes.

Out: chembl_activities.csv  ->  compound_id, smiles, gene, target_id,
                                target_name, act_type, pchembl

Run: conda activate npbs && python 12_get_chembl_api.py
"""
import csv, sys, time, os
import requests

HERE = os.path.dirname(os.path.abspath(__file__))

API = "https://www.ebi.ac.uk/chembl/api/data"
S = requests.Session()
S.headers["Accept"] = "application/json"

# Pain + neuroinflammation + diabetic/chronic-disease targets. Edit freely.
PANEL = [
    # nociception / ion channels
    "TRPV1", "TRPA1", "SCN9A", "SCN10A", "CACNA2D1", "P2RX3", "NTRK1",
    # neuroinflammation
    "TNF", "IL1B", "IL6", "NLRP3", "PTGS2", "NOS2", "TLR4", "NFKB1",
    # analgesic GPCRs / receptors
    "OPRM1", "CNR2", "ADRB2", "ADORA1", "HTR1A", "GRM5",
    # diabetic neuropathy / metabolic
    "AKR1B1", "PPARG", "DPP4", "PTPN1", "GSK3B", "SIRT1", "MAPK14", "RAGE",
]


def get(url, **params):
    for _ in range(5):
        try:
            r = S.get(url, params=params, timeout=60)
            r.raise_for_status()
            return r.json()
        except Exception as e:
            print("   retry:", e); time.sleep(5)
    raise RuntimeError(f"failed: {url} {params}")


def n_activities(tid):
    """count of pchembl>=6 activities for a target (0 for empty stub entries)."""
    j = get(f"{API}/activity.json", target_chembl_id=tid, pchembl_value__gte=6, limit=1)
    return j["page_meta"]["total_count"]


def resolve(gene):
    """gene symbol -> human SINGLE PROTEIN target with the MOST bioactivity.
    (Picking the first hit lands on empty curated stubs, e.g. TNF/GSK3B.)"""
    j = get(f"{API}/target/search.json", q=gene, limit=25)
    cands = [t for t in j.get("targets", [])
             if t.get("organism") == "Homo sapiens"
             and t.get("target_type") == "SINGLE PROTEIN"]
    if not cands:
        return None, None
    # ponytail: pick busiest target; if a gene has two legit isoform targets this takes the busier one
    best = max(cands, key=lambda t: n_activities(t["target_chembl_id"]))
    return best["target_chembl_id"], best["pref_name"]


def activities(tid):
    """all pchembl>=6, exact-relation activities for a target (paginated)."""
    url = f"{API}/activity.json"
    # keep ALL rows that carry a pchembl value (actives AND inactives);
    # thresholds are applied later so we get real negatives, not assumed ones.
    params = dict(target_chembl_id=tid, pchembl_value__isnull="false",
                  standard_relation="=", limit=1000)
    nxt = None
    while True:
        j = get(url if nxt is None else "https://www.ebi.ac.uk" + nxt, **({} if nxt else params))
        for a in j["activities"]:
            yield a
        nxt = j["page_meta"]["next"]
        if not nxt:
            break


def run():
    n = 0
    with open(os.path.join(HERE, "chembl_activities.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["compound_id", "smiles", "gene", "target_id",
                    "target_name", "act_type", "pchembl"])
        for gene in PANEL:
            tid, name = resolve(gene)
            if not tid:
                print(f"  {gene}: no target found -- skipped"); continue
            c = 0
            for a in activities(tid):
                smi = a.get("canonical_smiles")
                if not smi:
                    continue
                w.writerow([a["molecule_chembl_id"], smi, gene, tid, name,
                            a.get("standard_type"), a.get("pchembl_value")])
                c += 1; n += 1
            print(f"  {gene:10s} {tid:14s} {c:6d} activities")
    print(f"done: {n} rows -> chembl_activities.csv")
    assert n > 0, "no activities pulled -- check network / panel"
    print("check ok")


if __name__ == "__main__":
    run()
