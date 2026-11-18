"""P3 extension (amendment 06, Part A): the walk down the try order of each class.

Reads public/screen/p3x.order.tsv, the state of each candidate (public/<run>/P3_DONE, P3_FAILED, <run>.lengths.tsv.gz)
and public/screen/p3x_manual_exclusions.tsv (run_accession, reason: the kit rule, applied by hand after the fetch).
Per class, down the try order:
  - stop at MAX_ACCEPT accepted runs or MAX_TRIED candidates tried;
  - a candidate of a project with CAP accepted runs is passed over (not tried);
  - a fetched candidate is accepted if its subsample has >= 1 Gbp and N50 >= 9,000 and it is not excluded by hand;
    otherwise it is rejected with the reason; a failed fetch is rejected;
  - a candidate not yet fetched is "pending". Pending candidates are treated as accepted for the cap and the stop, so a
    batch never contains a run that the one-by-one walk would not reach.
Writes public/screen/p3x_status.tsv (every candidate reached) and public/screen/p3x_final.tsv (accepted runs; columns as
the frozen P3_03_confirm.py reads them). Prints the pending runs as "RUN READS BASES" lines (the next batch).
usage: python -I P3X_02_walk.py [--screen DIR] [--public DIR]
"""
import argparse
import os

import numpy as np
import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
MAX_ACCEPT, MAX_TRIED, CAP = 6, 12, 2
MIN_BASES, MIN_N50 = 1e9, 9000


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
    order = pd.read_csv(f"{a.screen}/p3x.order.tsv", sep="\t", dtype={"read_count": "int64", "base_count": "int64"})
    ex = {}
    if os.path.exists(f"{a.screen}/p3x_manual_exclusions.tsv"):
        t = pd.read_csv(f"{a.screen}/p3x_manual_exclusions.tsv", sep="\t", dtype=str)
        ex = dict(zip(t.run_accession, t.reason))
    rows, batch = [], []
    for cls, g in order.groupby("class", sort=True):
        acc, tried, per_project = 0, 0, {}
        for r in g.sort_values("try_order").itertuples():
            if acc >= MAX_ACCEPT or tried >= MAX_TRIED:
                break
            if per_project.get(r.project, 0) >= CAP:
                continue
            st, n, b = state(r.run_accession, a.public)
            tried += 1
            if st == "pending":
                status = "pending"
                batch.append((r.run_accession, r.read_count, r.base_count))
            elif st == "failed":
                status = "rejected: fetch failed"
            elif r.run_accession in ex:
                status = f"rejected: {ex[r.run_accession]}"
            elif b < MIN_BASES:
                status = "rejected: subsample < 1 Gbp"
            elif n < MIN_N50:
                status = "rejected: subsample N50 < 9 kb"
            else:
                status = "accepted"
            if status in ("accepted", "pending"):
                acc += 1
                per_project[r.project] = per_project.get(r.project, 0) + 1
            rows.append(dict(run_accession=r.run_accession, env=cls, project=r.project, unit=r.unit, try_order=r.try_order,
                             archive_mean_len=round(r.mean_len), status=status, subsample_n50=n, subsample_bases=b))
    s = pd.DataFrame(rows)
    s.to_csv(f"{a.screen}/p3x_status.tsv", sep="\t", index=False)
    s[s.status == "accepted"][["run_accession", "env", "project", "unit", "try_order", "subsample_n50", "subsample_bases"]].to_csv(
        f"{a.screen}/p3x_final.tsv", sep="\t", index=False)
    for run, n, b in batch:
        print(run, n, b)


if __name__ == "__main__":
    main()
