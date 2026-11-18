"""P3 step 3.2: status of the confirmatory runs and the replacement rule of amendment 05 Part A.

Run list public/screen/p3_runlist.tsv (created from p3_selected.tsv on the first call): run_accession, env, project,
unit, read_count, base_count, role ("selected" or "replaces <run>").
Status per run: ok | pending | fetch_failed (public/<run>/P3_FAILED) | short (subsample < 1 Gbp) | low_n50 (< 2 kb).
A run with a failing status that has no replacement yet gets the next unused unit of its project's replacement list
(public/screen/p3_replacements.tsv, in its seeded order). The kit rule is checked by hand and recorded in
public/screen/p3_manual_exclusions.tsv (run_accession, reason), which this script treats like a failing status.
Writes public/screen/p3_status.tsv and, for the analysis, public/screen/p3_final.tsv (the ok runs). Prints the number of
runs added.
usage: python -I P3_02b_check.py
"""
import os

import numpy as np
import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
S = os.environ.get("P3_SCREEN", f"{L}/public/screen")   # the override is for tests only
COLS = ["run_accession", "env", "project", "unit", "read_count", "base_count", "role"]


def n50(x):
    x = np.sort(np.asarray(x))[::-1]
    c = np.cumsum(x)
    return int(x[np.searchsorted(c, c[-1] / 2)])


def main():
    if not os.path.exists(f"{S}/p3_runlist.tsv"):
        pd.read_csv(f"{S}/p3_selected.tsv", sep="\t").assign(role="selected")[COLS].to_csv(
            f"{S}/p3_runlist.tsv", sep="\t", index=False)
    rl = pd.read_csv(f"{S}/p3_runlist.tsv", sep="\t")
    rep = pd.read_csv(f"{S}/p3_replacements.tsv", sep="\t")
    manual = {}
    if os.path.exists(f"{S}/p3_manual_exclusions.tsv"):
        m = pd.read_csv(f"{S}/p3_manual_exclusions.tsv", sep="\t")
        manual = dict(zip(m.run_accession, m.reason))
    rows = []
    for r in rl.itertuples():
        P = f"{L}/public/{r.run_accession}"
        st, nn, bases = "pending", np.nan, np.nan
        if r.run_accession in manual:
            st = f"excluded: {manual[r.run_accession]}"
        elif os.path.exists(f"{P}/P3_FAILED"):
            st = "fetch_failed"
        elif os.path.exists(f"{P}/P3_DONE"):
            ln = pd.read_csv(f"{P}/{r.run_accession}.lengths.tsv.gz", sep="\t").length
            nn, bases = n50(ln), int(ln.sum())
            st = "short" if bases < 1e9 else "low_n50" if nn < 2000 else "ok"
        rows.append(dict(run_accession=r.run_accession, status=st, subsample_n50=nn, subsample_bases=bases))
    stt = rl.merge(pd.DataFrame(rows), on="run_accession")
    replaced = set(stt.role[stt.role.str.startswith("replaces")].str.split().str[1])
    used_units = set(stt.unit)
    added = []
    for r in stt[~stt.status.isin(["ok", "pending"]) & ~stt.run_accession.isin(replaced)].itertuples():
        cand = rep[(rep.project == r.project) & ~rep.unit.isin(used_units)].sort_values("order")
        if cand.empty:
            continue
        c = cand.iloc[0]
        used_units.add(c.unit)
        added.append(dict(run_accession=c.run_accession, env=c.env, project=c.project, unit=c.unit,
                          read_count=int(c.read_count), base_count=int(c.base_count),
                          role=f"replaces {r.run_accession} ({r.status})"))
    if added:
        rl = pd.concat([rl, pd.DataFrame(added)[COLS]], ignore_index=True)
        rl.to_csv(f"{S}/p3_runlist.tsv", sep="\t", index=False)
    stt.to_csv(f"{S}/p3_status.tsv", sep="\t", index=False)
    stt[stt.status == "ok"][["run_accession", "env", "project", "unit", "role", "subsample_n50",
                             "subsample_bases"]].to_csv(f"{S}/p3_final.tsv", sep="\t", index=False)
    print(len(added))


if __name__ == "__main__":
    main()
