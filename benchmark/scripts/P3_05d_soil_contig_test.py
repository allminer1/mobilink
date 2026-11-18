"""P3 step 3.4a: the contig-level test on the deep-soil assemblies (per soil; no site comparison, vko8soil rule 8).

Molecule = contig; classes as registered (contigs/<assembly>.core.features.tsv.gz: ARG and MRG primary rule, MGE =
mobility core), annotated with the read pipeline on windows (P3_05c). Host stratum = the Binette bin of the contig
(vko8soil/bins/<assembly>/binette/final_contig_to_bin.tsv); other contigs are "unresolved". All A units of a contig
are one clump (P3_05b contig_clumps; calibrated in the no-linkage simulation P3_05e). Test: linkstat.pairs.net_test_pairs as amendment 04, at w = 2 kb (registered)
and 10 kb, with and without host strata.
Class sets: "registered", and "relaxed60" (ARG and MRG at identity >= 60%, the relaxed table of P0_02c; MGE unchanged),
because soil genes are far from the reference proteins and the registered classes leave few calls on contigs.
Next to it, per pair, the published pilot's statistic (contigs that carry both classes) under its label-permutation
null, on loci and on units (functions of P3_05b). Our soils are outside the registration.
usage: python -I P3_05d_soil_contig_test.py OUT_PREFIX [--assemblies mylo,flye] [--procs 8]
"""
import argparse
import importlib.util
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
V = "/mnt/rmrfc/workflow/iskander/vko8soil"
SOILS = ["278", "279", "280", "281", "304", "307", "308", "310"]
PAIRS = [("ARG", "MGE"), ("ARG", "MRG"), ("MRG", "MGE")]


def one(args):
    a, cs = args
    from linkstat.pairs import net_test_pairs

    spec = importlib.util.spec_from_file_location("p05b", f"{L}/scripts/P3_05b_pilot_test.py")
    pb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pb)
    ln = pd.read_csv(f"{L}/contigs/{a}.contigs.tsv", sep="\t").rename(columns={"contig": "molecule_id"})
    f = pd.read_csv(f"{L}/contigs/{a}.core.features.tsv.gz", sep="\t", dtype={"read_id": str}).rename(
        columns={"read_id": "molecule_id"})
    if cs == "relaxed60":
        r = pd.read_csv(f"{L}/contigs/{a}.relaxed.features.tsv.gz", sep="\t", dtype={"read_id": str}).rename(
            columns={"read_id": "molecule_id"})
        f = pd.concat([r[r["class"].isin(["ARG", "MRG"])], f[f["class"] == "MGE"]], ignore_index=True)
    f = f[f["class"].isin(["ARG", "MRG", "MGE"])][["molecule_id", "start", "end", "class", "scov"]].reset_index(drop=True)
    c2b = pd.read_csv(f"{V}/bins/{a}/binette/final_contig_to_bin.tsv", sep="\t", header=None,
                      names=["contig", "bin"]).set_index("contig").bin
    ln["host"] = ln.molecule_id.map(c2b).fillna("unresolved")
    f["host"] = f.molecule_id.map(ln.set_index("molecule_id").host)
    res, nul = [], []
    for w in (2000, 10000):
        for sn, st in (("bin", ln.host.to_numpy()), ("none", None)):
            t = net_test_pairs(ln[["molecule_id", "length"]], f, PAIRS, strata=st, w=w, oe=(w == 2000),
                               loci=pb.contig_clumps)
            t.insert(0, "strata", sn)
            t.insert(0, "class_set", cs)
            t.insert(0, "assembly", a)
            res.append(t)
    rng = np.random.default_rng(pb.SEED)
    for ca, cb in PAIRS:
        g = f[f["class"].isin([ca, cb])]
        for name, slots in (("loci", g), ("units", pb.to_units(g))):
            both = slots.groupby("molecule_id")["class"].nunique()
            obs = int((both > 1).sum())
            n_a, n_b = int((slots["class"] == ca).sum()), int((slots["class"] == cb).sum())
            nulls = {"C label permutation": pb.null_c(slots, rng, a=ca)}
            if name == "units":
                nulls["C label permutation within host"] = pb.null_c(slots, rng, by_host=True, a=ca)
            for nn, v in nulls.items():
                nul.append(dict(pb.summary(a, name, nn, obs, v, n_a, n_b), class_set=cs, class_a=ca, class_b=cb))
    base = dict(assembly=a, class_set=cs, contigs=len(ln), Mb=ln.length.sum() / 1e6, bins=int(c2b.nunique()),
                Mb_binned=ln.length[ln.host != "unresolved"].sum() / 1e6,
                **{f"calls_{c}": int((f["class"] == c).sum()) for c in ("ARG", "MRG", "MGE")},
                **{f"calls_{c}_binned": int(((f["class"] == c) & (f.host != "unresolved")).sum()) for c in ("ARG", "MRG", "MGE")})
    print(f"[3.4a] {a} {cs} done", flush=True)
    return pd.concat(res, ignore_index=True), pd.DataFrame(nul), base


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--assemblies", default="mylo,flye")
    ap.add_argument("--procs", type=int, default=8)
    x = ap.parse_args()
    todo = [(f"{p}_{s}", cs) for p in x.assemblies.split(",") for s in SOILS for cs in ("registered", "relaxed60")
            if os.path.exists(f"{L}/contigs/{p}_{s}.core.features.tsv.gz")]
    with Pool(x.procs) as p:
        r = p.map(one, todo, chunksize=1)
    pd.concat([a for a, _, _ in r], ignore_index=True).to_csv(f"{x.out}.linkstat.tsv", sep="\t", index=False)
    pd.concat([b for _, b, _ in r], ignore_index=True).rename(columns={"sample": "assembly"}).to_csv(
        f"{x.out}.nulls.tsv", sep="\t", index=False)
    pd.DataFrame([c for _, _, c in r]).to_csv(f"{x.out}.inputs.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
