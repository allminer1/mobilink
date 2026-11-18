"""P3 step 3.3: the sensitivity analyses of amendment 05 Part E (none changes a decision of Part D).

Each variant is the pooled registered test per class (as Part C: linkstat.pairs.pool_runs + net_test_pairs, no O/E),
with one thing changed. Loading follows the frozen P3_03_confirm.py (imported from frozen/amendment05).
Variants:
  primary            as P3_03 (a check: must equal tables/P3/confirm/pooled.tsv)
  chimera_filtered   reads flagged by linkstat.chimeras (overlap breakpoints, defaults) or with an adapter or barcode
                     inside the read (P2_05e table of P3_02, region "interior") removed
  genus_strata       Kraken2 strata at genus rank
  d500, d2000        unit distance
  ARG_MRG_60         ARG and MRG at identity >= 60% (the "relaxed" feature table)
  MGE_80             mobility-core MGE at identity >= 80%
  MGE_D3             MGE = the D3 classes of the primary feature table
  ARG_noint          ARG without the intrinsic templates
  no_few_hosts       without the runs whose per-run ARG:MGE result has fewer than 3 effective hosts
  R10_only, R9_only  by chemistry of the project (public/screen/p3_project_kits.tsv)
  loo:<run>          leave one run out
  with_PCR_runs, PCR_runs_only, pcr_run:<run>
                     not in Part E; added with deviation 1 (D66). The runs of p3_manual_exclusions.tsv that were
                     dropped as PCR-amplified libraries (PRJEB108131, human gut): pooled with the registered runs
                     of their class, pooled alone, and each run alone. Reported apart from the registered result.
usage: python -I P3_04_sensitivity.py RUN_LIST PER_RUN_TSV OUT_TSV [--procs 8] [--variants a,b,...]
"""
import argparse
import importlib.util
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
BASE = ["primary", "chimera_filtered", "genus_strata", "d500", "d2000", "ARG_MRG_60", "MGE_80", "MGE_D3", "ARG_noint",
        "no_few_hosts", "R10_only", "R9_only"]


def frozen():
    spec = importlib.util.spec_from_file_location("p303", f"{L}/frozen/amendment05/P3_03_confirm.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def load_variant(run, variant):
    """(lengths, features, strata, paf, exclude) for one run under a variant."""
    p303 = frozen()
    from linkstat.hosts import read_taxonomy, strata_for_molecules
    from linkstat.loci import read_paf

    ln, f, st, paf = p303.load_run(run)
    P = f"{L}/public/{run}/{run}"
    exclude = None
    if variant == "genus_strata":
        k = pd.read_csv(f"{P}.kraken.tsv.gz", sep="\t", dtype={"read_id": str})
        st = strata_for_molecules(ln.molecule_id.to_numpy(), k, read_taxonomy(f"{p303.KDB}/ktaxonomy.tsv"), rank="G")
    elif variant in ("ARG_MRG_60", "MGE_D3", "ARG_noint", "MGE_80"):
        mge = f[f["class"] == "MGE"]
        rest = f[f["class"] != "MGE"]
        if variant == "ARG_MRG_60":
            r = pd.read_csv(f"{P}.relaxed.features.tsv.gz", sep="\t", dtype={"read_id": str}).rename(
                columns={"read_id": "molecule_id"})
            rest = r[r["class"].isin(["ARG", "MRG"])][["molecule_id", "start", "end", "class", "scov"]]
        elif variant == "ARG_noint":
            r = pd.read_csv(f"{P}.primary.features.tsv.gz", sep="\t", dtype={"read_id": str}).rename(
                columns={"read_id": "molecule_id"})
            r = r[(r["class"] == "MRG") | ((r["class"] == "ARG") & ~r.intrinsic.astype(bool))]
            rest = r[["molecule_id", "start", "end", "class", "scov"]]
        elif variant == "MGE_D3":
            r = pd.read_csv(f"{P}.primary.features.tsv.gz", sep="\t", dtype={"read_id": str}).rename(
                columns={"read_id": "molecule_id"})
            mge = r[r["class"] == "MGE"][["molecule_id", "start", "end", "class", "scov"]]
        elif variant == "MGE_80":
            r = pd.read_csv(f"{P}.core.features.tsv.gz", sep="\t", dtype={"read_id": str}).rename(
                columns={"read_id": "molecule_id"})
            mge = r[(r["class"] == "MGE") & (r.pident >= 80)][["molecule_id", "start", "end", "class", "scov"]]
        f = pd.concat([rest, mge], ignore_index=True)
        f = f[f.molecule_id.isin(ln.molecule_id)].reset_index(drop=True)
        keep = set(f.loc[f["class"].isin(p303.A_CLASSES), "molecule_id"])
        paf = read_paf(f"{P}.ava.paf.gz", min_block=1000, keep=keep)
    elif variant == "chimera_filtered":
        from linkstat.chimeras import chimera_flags
        full = read_paf(f"{P}.ava.paf.gz", min_block=1000)
        fl = chimera_flags(paf=full)
        exclude = set(fl.loc[fl.flag, "read_id"])
        ad = f"{P}.adapters.tsv.gz"
        if os.path.exists(ad):
            a = pd.read_csv(ad, sep="\t", dtype={"read_id": str})
            exclude |= set(a.loc[a.region == "interior", "read_id"])
    return ln, f, st, paf, exclude


def job(args):
    env, runs, variant, out_dir = args
    path = f"{out_dir}/{env.replace(' ', '_')}__{variant.replace(':', '_')}.tsv"
    if os.path.exists(path):
        return path
    from linkstat.pairs import net_test_pairs, pool_runs
    p303 = frozen()
    base = variant.split(":")[0] if variant.startswith(("loo:", "pcr_run:")) else variant
    lv = "primary" if base in ("primary", "loo", "d500", "d2000", "no_few_hosts", "R10_only", "R9_only",
                              "with_PCR_runs", "PCR_runs_only", "pcr_run") else base
    data, n_excl = [], 0
    for run in runs:
        ln, f, st, paf, excl = load_variant(run, lv)
        if excl:
            keep = ~ln.molecule_id.isin(excl)
            n_excl += int((~keep).sum())
            st = np.asarray(st, dtype=object)[keep.to_numpy()]
            ln = ln[keep].reset_index(drop=True)
            f = f[f.molecule_id.isin(ln.molecule_id)]
            if paf is not None and len(paf):
                paf = paf[~paf.q.isin(excl) & ~paf.t.isin(excl)]
        data.append(dict(run=run, lengths=ln, features=f, strata=st, paf=paf))
    ln, f, st, paf = pool_runs(data)
    d = {"d500": 500, "d2000": 2000}.get(base, 1000)
    t = net_test_pairs(ln, f, p303.PAIRS, strata=st, paf=paf, oe=False, units_d=d)
    t.insert(0, "variant", variant)
    t.insert(0, "env", env)
    t["runs"], t["reads_removed"] = len(runs), n_excl
    t.to_csv(path + ".tmp", sep="\t", index=False)
    os.replace(path + ".tmp", path)
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_list")
    ap.add_argument("per_run")
    ap.add_argument("out")
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--variants", default=None)
    a = ap.parse_args()
    sel = pd.read_csv(a.run_list, sep="\t", dtype=str)
    pr = pd.read_csv(a.per_run, sep="\t")
    kits = pd.read_csv(f"{L}/public/screen/p3_project_kits.tsv", sep="\t", dtype=str)
    r9 = set(kits.loc[kits.chemistry.fillna("").str.startswith("R9"), "project"])
    few = set(pr.loc[(pr.class_a == "ARG") & (pr.class_b == "MGE") & (pr.hosts_effective < 3), "run"])
    out_dir = a.out + ".parts"
    os.makedirs(out_dir, exist_ok=True)
    variants = a.variants.split(",") if a.variants else BASE + ["loo", "pcr"]
    jobs = []
    ex = pd.read_csv(f"{L}/public/screen/p3_manual_exclusions.tsv", sep="\t", dtype=str)
    status = pd.read_csv(f"{L}/public/screen/p3_status.tsv", sep="\t", dtype=str).set_index("run_accession")
    pcr = [r for r in ex.loc[ex.reason.str.startswith("PCR-amplified"), "run_accession"]
           if os.path.exists(f"{L}/public/{r}/P3_DONE")]
    for env, g in sel.groupby("env"):
        runs = g.run_accession.tolist()
        for v in variants:
            if v == "pcr":
                mine = [r for r in pcr if status.env.get(r) == env]
                if mine:
                    jobs.append((env, runs + mine, "with_PCR_runs", out_dir))
                    jobs.append((env, mine, "PCR_runs_only", out_dir))
                    jobs += [(env, [r], f"pcr_run:{r}", out_dir) for r in mine]
            elif v == "loo":
                jobs += [(env, [x for x in runs if x != r], f"loo:{r}", out_dir) for r in runs if len(runs) > 1]
            elif v == "no_few_hosts":
                keep = [x for x in runs if x not in few]
                if keep and len(keep) < len(runs):
                    jobs.append((env, keep, v, out_dir))
            elif v in ("R10_only", "R9_only"):
                is9 = g.project.isin(r9)
                keep = g.run_accession[is9 if v == "R9_only" else ~is9].tolist()
                if keep and len(keep) < len(runs):
                    jobs.append((env, keep, v, out_dir))
            else:
                jobs.append((env, runs, v, out_dir))
    with Pool(a.procs) as p:
        paths = p.map(job, jobs, chunksize=1)
    t = pd.concat([pd.read_csv(x, sep="\t") for x in paths], ignore_index=True)
    t.to_csv(a.out, sep="\t", index=False)
    show = t[(t.class_a == "ARG") & ~t.variant.str.startswith(("loo:", "pcr_run:"))]
    print(show[["env", "variant", "class_b", "runs", "net_w", "net_w_lo", "net_w_hi", "net_w_p", "net_mol", "net_mol_p",
                "mol_ratio", "hosts_effective"]].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
