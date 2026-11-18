"""Registered classes for P2 analyses: ARG and MRG from the primary rule, MGE = the mobility core (D17).

MGE (registered, D17 / amendment 02): one call per ORF over all mobileOG-db hits at identity >= 60% and subject
cover >= 80% (P0_02c collapse), kept if the call's subject is in db/mobility_core.tsv.gz (as scripts/P1_06_g1b.py).
MGE_D3 (secondary): the primary-rule MGE calls (D3 categories) relabelled.
ARG, MRG and BIOCIDE: the primary-rule calls as they are.

usage: python -I P2_05i_core_features.py PRIMARY_FEATURES MOBILEOG_ALL OUT [--lift REF_CONTIGS]
  --lift: the inputs are 70 kb reference windows (P2_05b); the core MGE calls are lifted to contigs with the
          window-core rule of P2_05c, and PRIMARY_FEATURES is the already lifted ref.primary.features.tsv.gz.
"""
import importlib.util
import sys

import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    prim, mog, out = sys.argv[1], sys.argv[2], sys.argv[3]
    lift = sys.argv[5] if len(sys.argv) > 5 and sys.argv[4] == "--lift" else None
    p02c = load("p02c", f"{L}/scripts/P0_02c_features.py")
    d = p02c.collapse(p02c.read_hits(mog, 60.0))
    core = set(pd.read_csv(f"{L}/db/mobility_core.tsv.gz", sep="\t", usecols=["sseqid"]).sseqid)
    d = d[d.sseqid.isin(core)].rename(columns={"qseqid": "read_id"})
    mge = pd.DataFrame({"read_id": d.read_id.astype(str), "start": d.start.astype(int), "end": d.end.astype(int),
                        "class": "MGE", "gene": d.sseqid, "pident": d.pident, "scov": d.scov, "bitscore": d.bitscore})
    if lift:
        lf = load("lift", f"{L}/scripts/P2_05c_lift_ref_features.py")
        clen = pd.read_csv(lift, sep="\t").set_index("contig")["length"].to_dict()
        sp = mge.read_id.str.rsplit("@", n=1, expand=True)
        mge["contig"], mge["offset"] = sp[0], sp[1].astype(int)
        mge["clen"] = mge.contig.map(clen)
        lo = (mge.offset + lf.MARGIN).where(mge.offset > 0, 0)
        last = mge.offset + lf.WIN >= mge.clen
        hi = (mge.offset + lf.WIN - lf.MARGIN).where(~last, mge.clen + 1)
        mge = mge[(mge.start + mge.offset >= lo) & (mge.start + mge.offset < hi)].copy()
        mge["start"] += mge.offset
        mge["end"] += mge.offset
        mge["read_id"] = mge.contig
        mge = mge.drop(columns=["contig", "offset", "clen"])
    f = pd.read_csv(prim, sep="\t", dtype={"read_id": str})
    keep = f[f["class"].isin(["ARG", "MRG", "BIOCIDE"])]
    d3 = f[f["class"] == "MGE"].assign(**{"class": "MGE_D3"})
    cols = ["read_id", "start", "end", "class", "gene", "pident", "scov", "bitscore"]
    res = pd.concat([keep[cols], d3[cols], mge[cols]], ignore_index=True).sort_values(["read_id", "start", "class"])
    res.to_csv(out, sep="\t", index=False)
    print(f"[2.5] {out}: " + ", ".join(f"{k} {v}" for k, v in res["class"].value_counts().items()))


if __name__ == "__main__":
    main()
