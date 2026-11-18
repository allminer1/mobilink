#!/usr/bin/env python3
"""P0 step 0.2c: per-read feature tables from DIAMOND blastx hits, one soil at a time.

Classes
  ARG  CARD 4.0.1 protein homolog model. `intrinsic` flags vanR_in_vanO_cl, rpoB2 and HelR (vko8soil Phase 1).
  MRG  BacMet2 EXP entries with is_mrg == TRUE (vko8soil db/bacmet/bacmet_exp_classes.tsv). Entries without a
       metal are labelled BIOCIDE and are not used in the pair tests.
  MGE  mobileOG-db beatrix-1.6. Primary categories (user decision D3, 2026-11-14): integration/excision, transfer,
       replication/recombination/repair, phage. The fifth category, stability/transfer/defense, is labelled MGE_std.

Rules (identity %, subject cover %)
  primary  ARG 80/80  MRG 80/80  MGE 60/80   (docs/preregistration_H1-H4.md)
  strict   ARG 80/80  MRG 80/80  MGE 80/80
  relaxed  ARG 60/80  MRG 60/80  MGE 60/80
  mgeperm  ARG 80/80  MRG 80/80  MGE 40/80

Per rule and database: keep the hits that pass, then collapse overlapping hits on one read greedily by bitscore
(a hit is dropped if it overlaps an accepted hit of the same database by more than 50% of the shorter one). The
surviving hit gives the label, so one ORF gives one call per database.
Coordinates: DIAMOND qstart/qend are 1-based and qstart > qend on the minus strand; output is 0-based, end exclusive.

usage: P0_02c_features.py <soil> [--arg FILE] [--mrg FILE] [--mge FILE] [--out-dir DIR]
"""
import argparse
import gzip
import sys

import numpy as np
import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
V = "/mnt/rmrfc/workflow/iskander/vko8soil"
COLS = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "slen", "qlen", "evalue",
        "bitscore"]
RULES = {"primary": {"ARG": 80, "MRG": 80, "MGE": 60},   # pre-registration, fixed 2026-11-14
         "strict": {"ARG": 80, "MRG": 80, "MGE": 80},
         "relaxed": {"ARG": 60, "MRG": 60, "MGE": 60},
         "mgeperm": {"ARG": 80, "MRG": 80, "MGE": 40}}
MIN_SCOV = 80.0
INTRINSIC = {"vanR_in_vanO_cl", "rpoB2", "HelR"}
MGE_PRIMARY = {"integration/excision", "transfer", "replication/recombination/repair", "phage"}
OUT_COLS = ["read_id", "start", "end", "strand", "class", "gene", "label", "evidence", "intrinsic", "pident", "scov",
            "bitscore", "qlen", "db"]


def read_hits(path, min_id):
    """All hits with identity >= min_id and subject cover >= 80%, read in chunks. path 'none' skips the class."""
    parts = []
    if path == "none":
        path = None
    chunks = [] if path is None else pd.read_csv(path, sep="\t", names=COLS, header=None, chunksize=5_000_000,
                                                  dtype={"qseqid": str, "sseqid": str})
    for ch in chunks:
        scov = 100.0 * ((ch.send - ch.sstart).abs() + 1) / ch.slen
        keep = (ch.pident >= min_id) & (scov >= MIN_SCOV)
        ch = ch.loc[keep].copy()
        ch["scov"] = scov[keep]
        parts.append(ch)
    d = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=COLS + ["scov"])
    d["start"] = np.minimum(d.qstart, d.qend) - 1
    d["end"] = np.maximum(d.qstart, d.qend)
    d["strand"] = np.where(d.qstart <= d.qend, "+", "-")
    return d


def collapse(d):
    """Greedy by bitscore within each read: drop a hit overlapping an accepted one by > 50% of the shorter.
    Ties are broken by e-value, then subject ID, so the result does not depend on the input order."""
    if d.empty:
        return d
    d = d.sort_values(["qseqid", "bitscore", "evalue", "sseqid"], ascending=[True, False, True, True],
                      kind="stable").reset_index(drop=True)
    rid = d.qseqid.to_numpy(); s = d.start.to_numpy(); e = d.end.to_numpy()
    keep = np.zeros(len(d), dtype=bool)
    acc_s, acc_e, cur = [], [], None
    for i in range(len(d)):
        if rid[i] != cur:
            cur, acc_s, acc_e = rid[i], [], []
        ok = True
        for a, b in zip(acc_s, acc_e):
            ov = min(e[i], b) - max(s[i], a)
            if ov > 0.5 * min(e[i] - s[i], b - a):
                ok = False
                break
        if ok:
            keep[i] = True
            acc_s.append(s[i]); acc_e.append(e[i])
    return d.loc[keep].reset_index(drop=True)


def label_arg(d):
    f = d.sseqid.str.split("|", expand=True)
    d["gene"] = f[3]; d["label"] = f[2]; d["evidence"] = "CARD_homolog"
    d["intrinsic"] = d.gene.isin(INTRINSIC); d["class"] = "ARG"; d["db"] = "card_homolog"
    return d


def bacmet_key(sseqid):
    """BacMet2 EXP reuses 5 IDs (BAC0291, BAC0358, BAC0705-0707) for two proteins each, so the key is ID + accession
    (field 4, or field 3 in headers with fewer than 5 fields), as in vko8soil scripts/01_step06_densities.R."""
    f = sseqid.str.split("|")
    acc = np.where(f.str.len() >= 5, f.str[3], f.str[2])
    return f.str[0], acc


def label_mrg(d, classes):
    d["label"], d["Accession"] = bacmet_key(d.sseqid)
    m = d.merge(classes, left_on=["label", "Accession"], right_on=["BacMet_ID", "Accession"], how="left",
                validate="many_to_one")
    if m.is_mrg.isna().any() or len(m) != len(d):
        sys.exit(f"{int(m.is_mrg.isna().sum())} BacMet hits have no class entry")
    m["gene"] = m.Gene_name; m["evidence"] = m.metals.fillna("")
    m["intrinsic"] = False; m["class"] = np.where(m.is_mrg.astype(bool), "MRG", "BIOCIDE"); m["db"] = "bacmet2_exp"
    return m


def label_mge(d, lookup):
    """Gene, category and evidence from db/mobileog_all.lookup.tsv.gz: DIAMOND cuts sseqid at the first space,
    which removes the evidence field from 462,024 of 775,257 beatrix-1.6 headers."""
    m = d.merge(lookup, on="sseqid", how="left")
    if m.category.isna().any():
        sys.exit(f"{int(m.category.isna().sum())} mobileOG hits are missing from the lookup")
    m["label"] = m.category
    m["intrinsic"] = False; m["class"] = np.where(m.label.isin(MGE_PRIMARY), "MGE", "MGE_std"); m["db"] = "mobileog_all"
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("soil")
    ap.add_argument("--arg"); ap.add_argument("--mrg"); ap.add_argument("--mge")
    ap.add_argument("--lengths"); ap.add_argument("--out-dir", default=f"{L}/features")
    a = ap.parse_args()
    s = a.soil
    src = {"ARG": a.arg or f"{L}/dmnd/{s}.card.all.tsv.gz",
           "MRG": a.mrg or f"{V}/dmnd/{s}.bacmet.all.tsv.gz",
           "MGE": a.mge or f"{L}/dmnd/{s}.mobileog.all.tsv.gz"}
    classes = pd.read_csv(f"{V}/db/bacmet/bacmet_exp_classes.tsv", sep="\t",
                          usecols=["BacMet_ID", "Accession", "Gene_name", "metals", "is_mrg"], keep_default_na=False)
    classes["is_mrg"] = classes.is_mrg.astype(str).str.upper().map({"TRUE": True, "FALSE": False})
    if classes.is_mrg.isna().any():
        sys.exit("bacmet_exp_classes.tsv: is_mrg is not TRUE/FALSE everywhere")
    lookup = pd.read_csv(f"{L}/db/mobileog_all.lookup.tsv.gz", sep="\t", dtype=str)

    lengths = pd.read_csv(a.lengths or f"{L}/lengths/{s}.lengths.tsv.gz", sep="\t", dtype={"read_id": str})
    lengths = pd.Series(lengths.length.values, index=lengths.read_id.values)

    raw = {}
    for cls, path in src.items():
        min_id = min(r[cls] for r in RULES.values())
        raw[cls] = read_hits(path, min_id)
        # every hit read must be in the lengths table with the same length
        ql = lengths.reindex(raw[cls].qseqid.values).to_numpy()
        bad_missing = int(np.isnan(ql).sum()); bad_len = int(np.sum(ql[~np.isnan(ql)] != raw[cls].qlen.values[~np.isnan(ql)]))
        if bad_missing or bad_len:
            sys.exit(f"{s} {cls}: {bad_missing} hit reads missing from lengths, {bad_len} with a different length")

    summary = []
    for rule, thr in RULES.items():
        out = []
        label = {"ARG": label_arg, "MRG": lambda x: label_mrg(x, classes), "MGE": lambda x: label_mge(x, lookup)}
        for cls in ("ARG", "MRG", "MGE"):
            d = collapse(raw[cls].loc[raw[cls].pident >= thr[cls]].copy())
            if d.empty:
                continue
            out.append(label[cls](d).rename(columns={"qseqid": "read_id"})[OUT_COLS])
        f = pd.concat(out, ignore_index=True)
        f.to_csv(f"{a.out_dir}/{s}.{rule}.features.tsv.gz", sep="\t", index=False, compression="gzip")
        for c, g in f.groupby("class"):
            summary.append(dict(soil=s, rule=rule, cls=c, features=len(g), reads=g.read_id.nunique(),
                                intrinsic=int(g.intrinsic.sum()), mean_len=round(float((g.end - g.start).mean()), 1)))
    pd.DataFrame(summary).to_csv(f"{a.out_dir}/{s}.summary.tsv", sep="\t", index=False)
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    main()
