"""P4 step 4.6, part b (descriptive): ARG:MGE net linked fraction within 2 kb by sub-class, per environment class.

Sub-classes (fixed from the call counts of P4_06a before any linkage statistic was computed; D112):
  MGE side (all registered ARG against one kind of MGE):
    mge_intexc, mge_transfer      major mobileOG category of the mobility core call
  ARG side (one kind of ARG against all registered MGE):
    mech_*                        CARD resistance mechanism: efflux, inactivation, alteration, protection, replacement;
                                  other = several mechanisms or reduced permeability
    fam_*                         CARD drug class merged to families (beta-lactam; MLS = macrolide, lincosamide,
                                  streptogramin; sulfone with sulfonamide); a gene with two or more families is
                                  multidrug; families with fewer than 5,000 calls in the 241 runs form fam_other
  all                             check only: must equal the registered class estimate
Estimate and delete-one-project replicates by the functions of the frozen script (frozen/amendment07/P4_03_confirm.py).
Nothing here is a registered test.

usage: python -I P4_06b_subclass_linkage.py OUT_DIR [--procs 16] [--subs a,b] [--envs "marine,soil"]   (OUT_DIR: tables/P4/sub)
"""
import argparse
import importlib.util
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy import stats

L = "/mnt/rmrfc/workflow/iskander/linkstat"
CONF = f"{L}/tables/P4/confirm"
CARD = "/mnt/rmrfc/workflow/iskander/kmaarg/db/card/aro_index.tsv"
COLS = ["molecule_id", "start", "end", "class", "scov"]
BL = {"carbapenem", "cephalosporin", "cephamycin", "monobactam", "penem", "penam", "penicillin beta-lactam"}
MLS = {"macrolide antibiotic", "lincosamide antibiotic", "streptogramin antibiotic", "streptogramin A antibiotic", "streptogramin B antibiotic"}
MECH = {"antibiotic efflux": "mech_efflux", "antibiotic inactivation": "mech_inactivation", "antibiotic target alteration": "mech_alteration",
        "antibiotic target protection": "mech_protection", "antibiotic target replacement": "mech_replacement"}
FAMS = ["multidrug", "tetracycline", "MLS", "aminoglycoside", "beta-lactam", "fluoroquinolone", "peptide"]
SUBS = (["mge_intexc", "mge_transfer"] + list(MECH.values()) + ["mech_other"]
        + ["fam_" + x.replace("-", "_") for x in FAMS] + ["fam_other"])


def frozen_p4():
    s = importlib.util.spec_from_file_location("p4_03_confirm", f"{L}/frozen/amendment07/P4_03_confirm.py")
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


def family(drug):
    out = set()
    for x in drug.split(";"):
        x = "sulfonamide antibiotic" if x == "sulfone antibiotic" else x
        out.add("beta-lactam" if x in BL else "MLS" if x in MLS else x.replace(" antibiotic", ""))
    f = next(iter(out)) if len(out) == 1 else "multidrug"
    return "fam_" + (f.replace("-", "_") if f in FAMS else "other")


def gene_table():
    c = pd.read_csv(CARD, sep="\t", dtype=str).rename(columns={"CARD Short Name": "gene", "Drug Class": "drug", "Resistance Mechanism": "mech"})
    c = c[["gene", "drug", "mech"]].drop_duplicates("gene")
    return pd.DataFrame({"gene": c.gene, "mech": c.mech.map(MECH).fillna("mech_other"), "fam": c.drug.map(family)}).set_index("gene")


def load_sub(fz, run, sub, genes):
    ln, _, st = fz.load_light(run)
    f = pd.read_csv(f"{L}/public/{run}/{run}.core.features.tsv.gz", sep="\t", dtype={"read_id": str}).rename(columns={"read_id": "molecule_id"})
    f = f[f["class"].isin(["ARG", "MRG", "MGE"]) & f.molecule_id.isin(ln.molecule_id)]
    arg, mge = f[f["class"] == "ARG"], f[f["class"] == "MGE"]
    if sub.startswith("mge_"):
        cat = {"mge_intexc": "integration/excision", "mge_transfer": "transfer"}[sub]
        mge = mge[mge.gene.str.split("|").str[3] == cat]
    elif sub != "all":
        lab = arg.gene.map(genes["mech" if sub.startswith("mech_") else "fam"])
        assert lab.notna().all(), run
        arg = arg[(lab == sub).to_numpy()]
    f = pd.concat([arg[COLS], mge[COLS]], ignore_index=True)
    return dict(run=run, lengths=ln, features=f, strata=st), len(arg), len(mge)


def class_job(args):
    env, sel, sub, out = args
    path = f"{out}/linkage/{sub}/{env.replace(' ', '_')}.tsv"
    if os.path.exists(path):
        return path
    fz, genes = frozen_p4(), gene_table()
    data, na, nb = {}, 0, 0
    for r in sel.run_accession:
        data[r], a, b = load_sub(fz, r, sub, genes)
        na, nb = na + a, nb + b
    proj = sel.groupby("project").run_accession.apply(list)

    def est(dd):
        if not dd or not sum((d["features"]["class"] == "ARG").sum() for d in dd):
            return {}
        try:
            return fz.raw_arg_mge(dd)
        except ValueError as e:  # every ARG call of the sub-class shares its ORF with an MGE call: no unit is left
            if "at least one array" not in str(e) and "no class A call" not in str(e):  # the second: package v1.0
                raise
            return dict(note="no ARG unit left after the same-ORF rule")

    rows = [dict(env=env, sub=sub, left_out="none", runs=len(data), calls_a=na, calls_b=nb, **est(list(data.values())))]
    for p, rr in proj.items():
        rows.append(dict(env=env, sub=sub, left_out=p, runs=len(data) - len(rr), **est([d for k, d in data.items() if k not in rr])))
    pd.DataFrame(rows).to_csv(path + ".tmp", sep="\t", index=False)
    os.replace(path + ".tmp", path)
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--procs", type=int, default=16)
    ap.add_argument("--subs", default=",".join(SUBS))
    ap.add_argument("--envs", default="")
    a = ap.parse_args()
    subs = a.subs.split(",")
    assert set(subs) <= set(SUBS + ["all"])
    sel = pd.read_csv(f"{CONF}/p4_final.tsv", sep="\t", dtype=str)[["run_accession", "env", "project"]]
    envs = a.envs.split(",") if a.envs else sorted(sel.env.unique())
    for s in subs:
        os.makedirs(f"{a.out}/linkage/{s}", exist_ok=True)
    size = sel.groupby("env").size()
    jobs = sorted([(e, sel[sel.env == e], s, a.out) for s in subs for e in envs], key=lambda j: -size[j[0]])
    with Pool(min(a.procs, len(jobs))) as p:
        for path in p.imap_unordered(class_job, jobs, chunksize=1):
            print("done", path, flush=True)

    fz = frozen_p4()
    base = pd.read_csv(f"{CONF}/jackknife.tsv", sep="\t")
    rows = []
    for s in subs:
        for e in envs:
            g = pd.read_csv(f"{a.out}/linkage/{s}/{e.replace(' ', '_')}.tsv", sep="\t")
            f0 = g[g.left_out == "none"].iloc[0]
            full = float(f0.get("net_w", np.nan))
            reps = g[g.left_out != "none"].net_w.to_numpy(float) if "net_w" in g else np.array([])
            se, k = fz.jack_se(full, reps)
            t = stats.t.ppf(0.975, k - 1) if k >= 2 else np.nan
            rows.append(dict(sub=s, env=e, runs=int(f0.runs), projects=len(g) - 1, calls_a=int(f0.calls_a), calls_b=int(f0.calls_b),
                             units=f0.get("units", np.nan), sides_w=f0.get("sides_w", np.nan), events_w=f0.get("events_w", np.nan),
                             net_w=full, se_jack=se, lo95=full - t * se, hi95=full + t * se, jackknife_reps=k,
                             registered_all=float(base[(base.env == e) & (base.left_out == "none")].net_w.iloc[0])))
    name = "subclass.classes.tsv" if set(subs) == set(SUBS) and not a.envs else "subclass.partial.tsv"
    pd.DataFrame(rows).to_csv(f"{a.out}/{name}", sep="\t", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 300)
    print(pd.DataFrame(rows).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
