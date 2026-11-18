#!/usr/bin/env python3
"""P1 step 1.8: apply the G1b rule across runs (amendment 02 A3) to the per-run table of the frozen Part B driver.

Inputs: --table (frozen/partB/P1_06_g1b.py output, mobility core), --selected (public/screen/g1b_selected.tsv, for the
environment of each run). G1b passes if all of:
  at least 6 runs count, from at least 3 environments;
  at least two-thirds of the counting runs meet G1b;
  the runs that meet G1b come from at least 3 environments.
If fewer than 6 runs count (or they come from fewer than 3 environments), G1b is INCONCLUSIVE.
Secondary: median relative change native -> 3 kb (and -> 2 kb) over counting runs with |native| >= 0.10; within +-10%
allows "approximately length-invariant above N50 3 kb".
"""
import argparse

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--selected", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    t = pd.read_csv(a.table, sep="\t").merge(pd.read_csv(a.selected, sep="\t")[["run_accession", "env"]],
                                             left_on="run", right_on="run_accession", how="left")
    c = t[t.counts]
    n_count, env_count = len(c), c.env.nunique()
    meet = c[c.meets]
    share = len(meet) / n_count if n_count else np.nan
    if n_count < 6 or env_count < 3:
        verdict = "INCONCLUSIVE"
    elif share >= 2 / 3 and meet.env.nunique() >= 3:
        verdict = "PASS"
    else:
        verdict = "FAIL"
    rel3 = c.rel_change_3kb.dropna(); rel2 = c.rel_change_2kb.dropna()
    summ = dict(runs=len(t), counting=n_count, counting_envs=env_count, meeting=len(meet),
                meeting_share=share, meeting_envs=meet.env.nunique(), verdict=verdict,
                median_rel_change_3kb=rel3.median() if len(rel3) else np.nan, n_rel=len(rel3),
                median_rel_change_2kb=rel2.median() if len(rel2) else np.nan,
                wording=("approximately length-invariant above N50 3 kb" if len(rel3) and abs(rel3.median()) <= 0.10
                         else "much less length-dependent than the naive share"))
    pd.DataFrame([summ]).to_csv(a.out, sep="\t", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(t[["run", "env", "intact_n50", "intact_net", "r3000_net", "diff", "tol", "naive_ratio", "counts", "stable",
                 "overlap", "meets"]].round(4).to_string(index=False))
    for k, v in summ.items():
        print(f"{k}: {v}")


if __name__ == "__main__":
    main()
