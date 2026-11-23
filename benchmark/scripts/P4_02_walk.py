"""P4 (amendment 07, Part A): the state of the sample during and after the data pass.

Reads public/screen/p4_units.tsv (all units, with `selected`), public/screen/p4_selected.tsv (the sample), the state of
each run (public/<run>/P3_DONE, P3_FAILED, <run>.lengths.tsv.gz; the frozen run script writes these markers) and
public/screen/p4_manual_exclusions.tsv (run_accession, reason: the kit rule, applied by hand after the fetch).
  - a fetched run is accepted if its subsample has >= 1 Gbp and a read N50 >= 2,000 and it is not excluded by hand;
    otherwise it is rejected with the reason; a failed fetch is rejected;
  - a rejected run is replaced by the unit of the same project and class that was not selected and not yet tried and
    whose rank by archive mean read length is nearest to the rejected one (ties: the longer reads); if there is none,
    nothing replaces it;
  - a run not yet fetched is "pending".
Writes public/screen/p4_status.tsv (every run reached) and public/screen/p4_final.tsv (accepted runs: run_accession, env,
project, unit, subsample_n50, subsample_bases). Prints the pending runs as "RUN READS BASES" lines (the next batch).
usage: python -I P4_02_walk.py [--screen DIR] [--public DIR]
"""
import argparse
import os

import numpy as np
import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
MIN_BASES, MIN_N50 = 1e9, 2000


def n50(x):
    x = np.sort(np.asarray(x))[::-1]
    c = np.cumsum(x)
    return int(x[np.searchsorted(c, c[-1] / 2)])


def state(run, pub):
    d = f"{pub}/{run}"
    if os.path.exists(f"{d}/P3_FAILED"):
        return "failed", np.nan, np.nan
    if os.path.exists(f"{d}/P3_DONE"):
        ln = pd.read_csv(f"{d}/{run}.lengths.tsv.gz", sep="\t", usecols=["length"]).length
        ln = ln[ln > 0]
        return "done", n50(ln), int(ln.sum())
    return "pending", np.nan, np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--screen", default=f"{L}/public/screen")
    ap.add_argument("--public", default=f"{L}/public")
    a = ap.parse_args()
    U = pd.read_csv(f"{a.screen}/p4_units.tsv", sep="\t", dtype=str)
    U["mean_len"] = pd.to_numeric(U.mean_len)
    sel = pd.read_csv(f"{a.screen}/p4_selected.tsv", sep="\t", dtype=str)
    ena = pd.read_csv(f"{a.screen}/ena_p4_runs_2026-11-18.tsv", sep="\t", dtype=str).drop_duplicates("run_accession").set_index("run_accession")
    ex = {}
    if os.path.exists(f"{a.screen}/p4_manual_exclusions.tsv"):
        t = pd.read_csv(f"{a.screen}/p4_manual_exclusions.tsv", sep="\t", dtype=str)
        ex = dict(zip(t.run_accession, t.reason))
    rows, batch = [], []
    for (p, c), g in U[U["class"].isin(sel["class"].unique())].groupby(["study_accession", "class"], sort=True):
        g = g.sort_values(["mean_len", "run_accession"], ascending=[False, True]).reset_index(drop=True)
        tried = set()
        queue = [(i, "selected") for i in g.index[g.selected == "yes"]]
        while queue:
            i, role = queue.pop(0)
            r = g.loc[i]
            tried.add(i)
            st, n, b = state(r.run_accession, a.public)
            if st == "pending":
                status = "pending"
                batch.append((r.run_accession, ena.read_count[r.run_accession], ena.base_count[r.run_accession]))
            elif st == "failed":
                status = "rejected: fetch failed"
            elif r.run_accession in ex:
                status = f"rejected: {ex[r.run_accession]}"
            elif b < MIN_BASES:
                status = f"rejected: subsample {b / 1e9:.2f} Gbp"
            elif n < MIN_N50:
                status = f"rejected: subsample N50 {n}"
            else:
                status = "accepted"
            rows.append(dict(run_accession=r.run_accession, env=c, project=p, unit=r.unit, rank=i + 1, role=role, status=status,
                             subsample_n50=n, subsample_bases=b))
            if status.startswith("rejected"):
                free = [j for j in g.index if j not in tried and g.selected[j] != "yes" and all(j != q for q, _ in queue)]
                if free:
                    j = min(free, key=lambda k: (abs(k - i), k))
                    queue.append((j, f"replaces {r.run_accession}"))
    S = pd.DataFrame(rows)
    S.to_csv(f"{a.screen}/p4_status.tsv", sep="\t", index=False)
    S[S.status == "accepted"][["run_accession", "env", "project", "unit", "subsample_n50", "subsample_bases"]].to_csv(
        f"{a.screen}/p4_final.tsv", sep="\t", index=False)
    for r in batch:
        print(*r)


if __name__ == "__main__":
    main()
