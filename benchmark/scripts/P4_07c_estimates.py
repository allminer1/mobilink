"""P4 steps 4.7 (Part 2) and 4.7b (exploratory; D114, D115): naive share and net linked fraction per run and per group
of runs, with the registered classes, for a table of runs annotated by the frozen run script.

Per run and per group (runs pooled as in amendment 05 Part C: anchors stay in their run and host stratum), for the
pairs ARG:MGE and MRG:MGE:
  naive     share of the A-carrying reads that also carry a B call (after the same-ORF rule);
  net_w     net linked fraction within 2 kb (mobilink.pairs.net_test_pairs without unit loci, 99 draws: estimate, SE and
            95% interval of the package); for groups of two or more runs also the delete-one-run jackknife SE of the
            estimate (functions of the frozen script of amendment 07).
Nothing here is a registered test.

usage: python -I P4_07c_estimates.py RUN_TABLE GROUP_COLUMNS OUT_PREFIX [--procs 8]
  RUN_TABLE      tab-separated, with run_accession and the group columns (e.g. public/screen/p4_07_runs.tsv)
  GROUP_COLUMNS  comma-separated (e.g. source   or   set,sample,condition)
"""
import argparse
import importlib.util
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
PAIRS = [("ARG", "MGE"), ("MRG", "MGE")]
KEEP = ["units", "sides_w", "net_w", "net_w_se", "net_w_lo", "net_w_hi", "net_w_p_enrich"]


def frozen_p4():
    s = importlib.util.spec_from_file_location("p4_03_confirm", f"{L}/frozen/amendment07/P4_03_confirm.py")
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


def naive(f, a, b):
    from mobilink.features import drop_same_orf

    fa, fb, _ = drop_same_orf(f[f["class"] == a], f[f["class"] == b], 0.5)
    ra = fa.molecule_id.unique()
    return len(ra), int(np.isin(ra, fb.molecule_id.unique()).sum())


def estimates(ln, f, st):
    from mobilink.fragment import n50
    from mobilink.pairs import net_test_pairs

    x = ln.length.to_numpy(np.int64)
    row = dict(reads=len(x), bases=int(x.sum()), n50=n50(x), mean_len=float(x.mean()), **{f"calls_{c}": int((f["class"] == c).sum()) for c in ("ARG", "MRG", "MGE")})
    for a, b in PAIRS:
        na, nab = naive(f, a, b)
        row.update({f"{a}_{b}_reads_A": na, f"{a}_{b}_naive": nab / na if na else np.nan})
        if na:
            try:
                t = net_test_pairs(ln, f, [(a, b)], strata=st, loci=None, n_mc=99, oe=False).iloc[0]
                row.update({f"{a}_{b}_{k}": t.get(k) for k in KEEP})
            except Exception as e:  # too few calls for the estimator
                row[f"{a}_{b}_note"] = str(e)[:80]
    return row


def group_job(args):
    key, runs, cols = args
    from mobilink.pairs import pool_runs

    fz = frozen_p4()
    data, out_runs = [], []
    for r in runs:
        if not os.path.exists(f"{L}/public/{r}/P3_DONE"):
            out_runs.append(dict(run_accession=r, fetched=False))
            continue
        ln, f, st = fz.load_light(r)
        data.append(dict(run=r, lengths=ln, features=f, strata=st))
        out_runs.append(dict(run_accession=r, fetched=True, **estimates(ln, f, st)))
    g = dict(zip(cols, key), runs_listed=len(runs), runs=len(data))
    if data:
        ln, f, st, _ = pool_runs(data)
        g.update(estimates(ln, f, st))
        if len(data) >= 2:
            for a, b in PAIRS:
                reps = []
                for i in range(len(data)):
                    l2, f2, s2, _ = pool_runs(data[:i] + data[i + 1:])
                    fa, fb = f2[f2["class"] == a], f2[f2["class"] == b]
                    reps.append(fz.fast_net_w(l2, fa, fb, s2)["net_w"] if len(fa) and len(fb) else np.nan)
                se, k = fz.jack_se(g.get(f"{a}_{b}_net_w", np.nan), reps)
                g.update({f"{a}_{b}_net_w_se_jack": se, f"{a}_{b}_jack_reps": k})
    return g, out_runs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("group")
    ap.add_argument("out")
    ap.add_argument("--procs", type=int, default=8)
    a = ap.parse_args()
    cols = a.group.split(",")
    t = pd.read_csv(a.table, sep="\t", dtype=str)
    jobs = [(k if isinstance(k, tuple) else (k,), g.run_accession.tolist(), cols) for k, g in t.groupby(cols if len(cols) > 1 else cols[0])]
    jobs.sort(key=lambda j: -len(j[1]))
    with Pool(min(a.procs, len(jobs))) as p:
        res = p.map(group_job, jobs, chunksize=1)
    pd.DataFrame([g for g, _ in res]).sort_values(cols).to_csv(f"{a.out}.groups.tsv", sep="\t", index=False)
    r = pd.DataFrame([x for _, rr in res for x in rr])
    t.merge(r, on="run_accession", how="left").to_csv(f"{a.out}.runs.tsv", sep="\t", index=False)
    print("groups", len(res), "| runs fetched", int(r.fetched.sum()), "of", len(r))


if __name__ == "__main__":
    main()
