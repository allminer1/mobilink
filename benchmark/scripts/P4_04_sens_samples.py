"""P4 step 4.5b (amendment 07 Part E and note 4): the sensitivity analyses that only change the sample.

Scenarios (ARG:MGE class estimates, delete-one-project jackknife and the two H3a contrasts, computed with the
functions of the frozen script frozen/amendment07/P4_03_confirm.py; no cut, no pooled test):
  pre2021    runs made public before 2021 left out (Part E);
  danish     the projects of the Danish consortium left out together (Part E): PRJEB104930, PRJEB114667, PRJEB85558
             (the fourth, PRJEB125553, has no accepted run);
  megabase   runs in which reads over 200 kb hold more than 1% of the bases left out (note 4; exploratory).
Classes that a scenario does not touch keep the jackknife table of the confirmatory run.
None of these changes a registered verdict; they are reported next to it.

usage: python -I P4_04_sens_samples.py OUT_DIR [--procs 6]     (OUT_DIR: tables/P4/sens)
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
DANISH = ["PRJEB104930", "PRJEB114667", "PRJEB85558", "PRJEB125553"]
H3A = [("wastewater", "soil"), ("human gut", "soil")]


def frozen_p4():
    s = importlib.util.spec_from_file_location("p4_03_confirm", f"{L}/frozen/amendment07/P4_03_confirm.py")
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


def job(args):
    env, sel, out = args
    return frozen_p4().class_jack((env, sel, out, True))


def megabase_runs(sel):
    out = []
    for r in sel.run_accession:
        x = pd.read_csv(f"{L}/public/{r}/{r}.lengths.tsv.gz", sep="\t", usecols=["length"]).length.to_numpy(np.int64)
        if x[x > 200000].sum() > 0.01 * x.sum():
            out.append(r)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--procs", type=int, default=6)
    a = ap.parse_args()
    sel = pd.read_csv(f"{CONF}/p4_final.tsv", sep="\t", dtype=str)[["run_accession", "env", "project"]]
    ena = pd.read_csv(f"{L}/public/screen/ena_p4_runs_2026-11-18.tsv", sep="\t", dtype=str, usecols=["run_accession", "first_public"])
    first = sel.merge(ena, on="run_accession", how="left").first_public
    assert first.notna().all()
    drop = {"pre2021": sel.run_accession[(first.str[:4] < "2021").to_numpy()].tolist(),
            "danish": sel.run_accession[sel.project.isin(DANISH)].tolist(),
            "megabase": megabase_runs(sel)}
    jobs, keep = [], {}
    for name, runs in drop.items():
        o = f"{a.out}/{name}"
        for d in ("jackknife", "projects", "cut"):
            os.makedirs(f"{o}/{d}", exist_ok=True)
        pd.DataFrame({"run_accession": runs}).merge(sel).to_csv(f"{o}/left_out.tsv", sep="\t", index=False)
        keep[name] = sel[~sel.run_accession.isin(runs)]
        for env in sorted(sel.env[sel.run_accession.isin(runs)].unique()):
            jobs.append((env, keep[name][keep[name].env == env], o))
        print(name, len(runs), "runs left out;", sorted(sel.env[sel.run_accession.isin(runs)].unique()), flush=True)
    with Pool(min(a.procs, len(jobs))) as p:
        p.map(job, jobs, chunksize=1)

    fz = frozen_p4()
    base = pd.read_csv(f"{CONF}/jackknife.tsv", sep="\t")
    rows, cons = [], []
    for name in ["registered"] + list(drop):
        s = sel if name == "registered" else keep[name]
        C = {}
        for env in sorted(s.env.unique()):
            path = f"{a.out}/{name}/jackknife/{env.replace(' ', '_')}.tsv"
            g = pd.read_csv(path, sep="\t") if os.path.exists(path) else base[base.env == env]
            full = float(g[g.left_out == "none"].net_w.iloc[0])
            se, k = fz.jack_se(full, g[g.left_out != "none"].net_w)
            t = stats.t.ppf(0.975, k - 1) if k >= 2 else np.nan
            C[env] = dict(net_w=full, se=se, df=k - 1, reps=g[g.left_out != "none"].net_w.to_numpy(float))
            rows.append(dict(scenario=name, env=env, recomputed=os.path.exists(path), runs=int((s.env == env).sum()),
                             projects=s[s.env == env].project.nunique(), units=int(g[g.left_out == "none"].units.iloc[0]),
                             net_w=full, se_jack=se, lo95=full - t * se, hi95=full + t * se))
        for x, y in H3A:
            ra, rb = C[x], C[y]
            d, se = ra["net_w"] - rb["net_w"], float(np.hypot(ra["se"], rb["se"]))
            df = se ** 4 / (ra["se"] ** 4 / ra["df"] + rb["se"] ** 4 / rb["df"])
            cons.append(dict(scenario=name, contrast=f"{x} - {y}", diff=d, se=se, df=df, t=d / se, p_one_sided=float(stats.t.sf(d / se, df)),
                             smallest_without_one_project=float(min(np.nanmin(ra["reps"]) - rb["net_w"], ra["net_w"] - np.nanmax(rb["reps"]))),
                             largest_without_one_project=float(max(np.nanmax(ra["reps"]) - rb["net_w"], ra["net_w"] - np.nanmin(rb["reps"])))))
    pd.DataFrame(rows).to_csv(f"{a.out}/sens_samples.classes.tsv", sep="\t", index=False)
    pd.DataFrame(cons).to_csv(f"{a.out}/sens_samples.h3a.tsv", sep="\t", index=False)
    print(pd.DataFrame(cons).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
