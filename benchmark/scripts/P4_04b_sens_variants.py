"""P4 step 4.5b (amendment 07 Part E): the sensitivity analyses that change the reads, the strata or the classes.

For each variant and class: the pooled ARG:MGE class estimate and its delete-one-project jackknife, computed as in
the frozen script (frozen/amendment07/P4_03_confirm.py: load_light, raw_arg_mge, jack_se), with one thing changed.
No cut, no pooled test. None of these changes a registered verdict; they are reported next to it.
Variants:
  primary            no change (a check: must equal tables/P4/confirm/jackknife.tsv)
  chimera_filtered   reads flagged by mobilink.chimeras (overlap breakpoints, defaults) or with an adapter or barcode
                     inside the read (adapter table of the data pass, region "interior") removed
  control_removed    suspected control-strand reads removed: Kraken2 taxid 562 and length 3,400-3,619, in the runs
                     that the screen of D92 flags (>= 100 such reads and >= 30% of the reads labelled 562), or where
                     >= 100 such reads are >= 1.5 times the mean of the two neighbouring length windows among the
                     reads labelled 562 (the second criterion as applied in P3, D92)
  genus_strata       Kraken2 strata at genus rank
  ARG_MRG_60         ARG and MRG at identity >= 60% (the "relaxed" feature table)
  MGE_80             mobility-core MGE at identity >= 80%
  MGE_D3             MGE = the broad D3 classes of the primary feature table
Output: OUT_DIR/variants/<variant>/<class>.tsv, OUT_DIR/sens_variants.classes.tsv, OUT_DIR/sens_variants.h3a.tsv,
OUT_DIR/control_screen.tsv (per run), OUT_DIR/chimera_counts.tsv (per run).
usage: python -I P4_04b_sens_variants.py OUT_DIR [--procs 12] [--variants a,b]     (OUT_DIR: tables/P4/sens)
"""
import argparse
import gzip
import importlib.util
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy import stats

L = "/mnt/rmrfc/workflow/iskander/linkstat"
CONF = f"{L}/tables/P4/confirm"
VARIANTS = ["primary", "chimera_filtered", "control_removed", "genus_strata", "ARG_MRG_60", "MGE_80", "MGE_D3"]
H3A = [("wastewater", "soil"), ("human gut", "soil")]
LO, HI = 3400, 3620
COLS = ["molecule_id", "start", "end", "class", "scov"]


def frozen_p4():
    s = importlib.util.spec_from_file_location("p4_03_confirm", f"{L}/frozen/amendment07/P4_03_confirm.py")
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


def excl_path(out, kind, run):
    return f"{out}/excl/{kind}/{run}.txt.gz"


def write_ids(path, ids):
    with gzip.open(path + ".tmp", "wt") as fh:
        fh.write("".join(f"{x}\n" for x in ids))
    os.replace(path + ".tmp", path)


def read_ids(path):
    with gzip.open(path, "rt") as fh:
        return {x.rstrip("\n") for x in fh if x.strip()}


def chimera_job(args):
    run, out = args
    path, P = excl_path(out, "chimera", run), f"{L}/public/{run}/{run}"
    if not os.path.exists(path):
        from mobilink.chimeras import chimera_flags
        from mobilink.loci import read_paf

        fl = chimera_flags(paf=read_paf(f"{P}.ava.paf.gz", min_block=1000))
        ids = set(fl.loc[fl.flag, "read_id"].astype(str))
        a = pd.read_csv(f"{P}.adapters.tsv.gz", sep="\t", dtype={"read_id": str})
        ad = set(a.loc[a.region == "interior", "read_id"])
        write_ids(path, sorted(ids | ad))
        pd.DataFrame([dict(run=run, judged=len(fl), flagged=len(ids), interior_adapter=len(ad), removed=len(ids | ad))]).to_csv(
            path.replace(".txt.gz", ".counts.tsv"), sep="\t", index=False)
    return pd.read_csv(path.replace(".txt.gz", ".counts.tsv"), sep="\t")


def control_job(args):
    run, out = args
    path, P = excl_path(out, "control", run), f"{L}/public/{run}/{run}"
    cpath = path.replace(".txt.gz", ".screen.tsv")
    if not os.path.exists(cpath):
        ln = pd.read_csv(f"{P}.lengths.tsv.gz", sep="\t", dtype={"read_id": str})
        k = pd.read_csv(f"{P}.kraken.tsv.gz", sep="\t", dtype={"read_id": str}, usecols=["read_id", "taxid"])
        e = ln[ln.read_id.isin(set(k.read_id[k.taxid == 562]))]
        x, w = e.length.to_numpy(), HI - LO
        win = (x >= LO) & (x < HI)
        nb = (((x >= LO - w) & (x < LO)).sum() + ((x >= HI) & (x < HI + w)).sum()) / 2
        ce = int(win.sum())
        ratio = float(ce / nb) if nb else np.nan
        share = ce / len(e) if len(e) else np.nan
        flag = bool(ce >= 100 and share >= 0.3)
        flag2 = bool(ce >= 100 and (np.isnan(ratio) or ratio >= 1.5))
        write_ids(path, sorted(e.read_id[win]) if (flag or flag2) else [])
        pd.DataFrame([dict(run=run, reads=len(ln), ecoli_reads=len(e), ecoli_window_reads=ce, ecoli_window_share=share,
                           ecoli_window_ratio=ratio, window_share_of_all_reads=ce / len(ln), flag_d92=flag, flag_ratio=flag2,
                           removed=ce if (flag or flag2) else 0)]).to_csv(cpath, sep="\t", index=False)
    return pd.read_csv(cpath, sep="\t")


def load_variant(fz, run, variant, out):
    from mobilink.hosts import read_taxonomy, strata_for_molecules

    ln, f, st = fz.load_light(run)
    P = f"{L}/public/{run}/{run}"
    if variant == "genus_strata":
        k = pd.read_csv(f"{P}.kraken.tsv.gz", sep="\t", dtype={"read_id": str})
        st = strata_for_molecules(ln.molecule_id.to_numpy(), k, read_taxonomy(f"{fz.KDB}/ktaxonomy.tsv"), rank="G")
    elif variant in ("ARG_MRG_60", "MGE_80", "MGE_D3"):
        rd = lambda name: pd.read_csv(f"{P}.{name}.features.tsv.gz", sep="\t", dtype={"read_id": str}).rename(columns={"read_id": "molecule_id"})  # noqa: E731
        mge, rest = f[f["class"] == "MGE"], f[f["class"] != "MGE"]
        if variant == "ARG_MRG_60":
            r = rd("relaxed")
            rest = r[r["class"].isin(["ARG", "MRG"])][COLS]
        elif variant == "MGE_D3":
            r = rd("primary")
            mge = r[r["class"] == "MGE"][COLS]
        else:
            r = rd("core")
            mge = r[(r["class"] == "MGE") & (r.pident >= 80)][COLS]
        f = pd.concat([rest[COLS], mge[COLS]], ignore_index=True)
        f = f[f.molecule_id.isin(ln.molecule_id)].reset_index(drop=True)
    elif variant in ("chimera_filtered", "control_removed"):
        ex = read_ids(excl_path(out, "chimera" if variant == "chimera_filtered" else "control", run))
        if ex:
            keep = ~ln.molecule_id.isin(ex).to_numpy()
            st = np.asarray(st, dtype=object)[keep]
            ln = ln[keep].reset_index(drop=True)
            f = f[f.molecule_id.isin(ln.molecule_id)].reset_index(drop=True)
    return dict(run=run, lengths=ln, features=f, strata=st)


def class_job(args):
    env, sel, variant, out = args
    path = f"{out}/variants/{variant}/{env.replace(' ', '_')}.tsv"
    if os.path.exists(path):
        return path
    fz = frozen_p4()
    data = {r: load_variant(fz, r, variant, out) for r in sel.run_accession}
    reads = int(sum(len(d["lengths"]) for d in data.values()))
    proj = sel.groupby("project").run_accession.apply(list)
    rows = [dict(env=env, variant=variant, left_out="none", runs=len(data), reads=reads, **fz.raw_arg_mge(list(data.values())))]
    for p, rr in proj.items():
        rest = [d for k, d in data.items() if k not in rr]
        rows.append(dict(env=env, variant=variant, left_out=p, runs=len(rest), **(fz.raw_arg_mge(rest) if rest else {})))
    pd.DataFrame(rows).to_csv(path + ".tmp", sep="\t", index=False)
    os.replace(path + ".tmp", path)
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--procs", type=int, default=12)
    ap.add_argument("--variants", default=",".join(VARIANTS))
    a = ap.parse_args()
    variants = a.variants.split(",")
    assert set(variants) <= set(VARIANTS)
    sel = pd.read_csv(f"{CONF}/p4_final.tsv", sep="\t", dtype=str)[["run_accession", "env", "project"]]
    for d in ["excl/chimera", "excl/control"] + [f"variants/{v}" for v in variants]:
        os.makedirs(f"{a.out}/{d}", exist_ok=True)
    order = sel.assign(s=[os.path.getsize(f"{L}/public/{r}/{r}.ava.paf.gz") for r in sel.run_accession]).sort_values("s", ascending=False)
    with Pool(a.procs) as p:
        if "control_removed" in variants:
            t = pd.concat(p.map(control_job, [(r, a.out) for r in sel.run_accession], chunksize=4), ignore_index=True).merge(sel, left_on="run", right_on="run_accession")
            t.drop(columns="run_accession").to_csv(f"{a.out}/control_screen.tsv", sep="\t", index=False)
            print("control screen: runs with removal", int((t.removed > 0).sum()), "| flag_d92", int(t.flag_d92.sum()), "| flag_ratio", int(t.flag_ratio.sum()), flush=True)
    if "chimera_filtered" in variants:
        with Pool(min(a.procs, 8)) as p:      # the overlap tables are held in memory
            t = pd.concat(p.map(chimera_job, [(r, a.out) for r in order.run_accession], chunksize=1), ignore_index=True).merge(sel, left_on="run", right_on="run_accession")
            t.drop(columns="run_accession").to_csv(f"{a.out}/chimera_counts.tsv", sep="\t", index=False)
            print("chimera flags done: median share of judged reads flagged", round(float((t.flagged / t.judged.clip(lower=1)).median()), 4), flush=True)
    envs = sorted(sel.env.unique())
    size = sel.groupby("env").size()
    jobs = sorted([(e, sel[sel.env == e], v, a.out) for v in variants for e in envs], key=lambda j: -size[j[0]])
    with Pool(a.procs) as p:
        p.map(class_job, jobs, chunksize=1)

    fz = frozen_p4()
    base = pd.read_csv(f"{CONF}/jackknife.tsv", sep="\t")
    rows, cons = [], []
    for v in variants:
        C = {}
        for e in envs:
            g = pd.read_csv(f"{a.out}/variants/{v}/{e.replace(' ', '_')}.tsv", sep="\t")
            full = float(g[g.left_out == "none"].net_w.iloc[0])
            reps = g[g.left_out != "none"].net_w.to_numpy(float)
            se, k = fz.jack_se(full, reps)
            t = stats.t.ppf(0.975, k - 1) if k >= 2 else np.nan
            reg = float(base[(base.env == e) & (base.left_out == "none")].net_w.iloc[0])
            C[e] = dict(net_w=full, se=se, df=k - 1, reps=reps)
            rows.append(dict(variant=v, env=e, runs=int(g[g.left_out == "none"].runs.iloc[0]), reads=g[g.left_out == "none"].reads.iloc[0],
                             units=g[g.left_out == "none"].units.iloc[0], net_w=full, se_jack=se, lo95=full - t * se, hi95=full + t * se,
                             jackknife_reps=k, registered=reg, change=full - reg))
        for x, y in H3A:
            ra, rb = C[x], C[y]
            if not (np.isfinite(ra["se"]) and np.isfinite(rb["se"])):
                cons.append(dict(variant=v, contrast=f"{x} - {y}"))
                continue
            d, se = ra["net_w"] - rb["net_w"], float(np.hypot(ra["se"], rb["se"]))
            df = se ** 4 / (ra["se"] ** 4 / ra["df"] + rb["se"] ** 4 / rb["df"])
            cons.append(dict(variant=v, contrast=f"{x} - {y}", diff=d, se=se, df=df, t=d / se, p_one_sided=float(stats.t.sf(d / se, df)),
                             smallest_without_one_project=float(min(np.nanmin(ra["reps"]) - rb["net_w"], ra["net_w"] - np.nanmax(rb["reps"]))),
                             largest_without_one_project=float(max(np.nanmax(ra["reps"]) - rb["net_w"], ra["net_w"] - np.nanmin(rb["reps"])))))
    pd.DataFrame(rows).to_csv(f"{a.out}/sens_variants.classes.tsv", sep="\t", index=False)
    pd.DataFrame(cons).to_csv(f"{a.out}/sens_variants.h3a.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(cons).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
