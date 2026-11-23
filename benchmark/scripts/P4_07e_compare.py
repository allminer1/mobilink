"""P4 step 4.7, Part 2 (exploratory; D114): the 40 runs of Go and Shrestha 2026 with our pipeline, per environment.

Inputs: OURS_PREFIX.groups.tsv / .runs.tsv from P4_07c_estimates.py (group column `source`), and the Part 1 tables
tables/P4/reanalysis/their_env.tsv and their_runs.tsv.
Outputs OUT_PREFIX.env.tsv (one row per environment: their share, our naive share, our net linked fraction within
2 kb with the package interval and the delete-one-run jackknife SE, read N50, amplified runs, and the rank by each
measure) and OUT_PREFIX.stats.tsv (rank correlations; agreement of their share with our naive share per run).
An interval from the delete-one-run jackknife uses t with (runs - 1) degrees of freedom. Runs of an environment come
from one or two projects, so these intervals describe the runs at hand and not the environment.

usage: python -I P4_07e_compare.py OURS_PREFIX OUT_PREFIX
"""
import sys

import numpy as np
import pandas as pd
from scipy import stats

L = "/mnt/rmrfc/workflow/iskander/linkstat"
R = f"{L}/tables/P4/reanalysis"


def main():
    ours, out = sys.argv[1], sys.argv[2]
    g = pd.read_csv(f"{ours}.groups.tsv", sep="\t")
    r = pd.read_csv(f"{ours}.runs.tsv", sep="\t")
    te = pd.read_csv(f"{R}/their_env.tsv", sep="\t")
    tr = pd.read_csv(f"{R}/their_runs.tsv", sep="\t")
    e = te[["source", "projects", "amplified_runs", "share_pooled", "cut_share_pooled"]].rename(
        columns={"share_pooled": "their_share", "cut_share_pooled": "their_share_cut_2kb"})
    e = e.merge(g[["source", "runs_listed", "runs", "reads", "n50", "calls_ARG", "calls_MGE", "ARG_MGE_reads_A", "ARG_MGE_naive", "ARG_MGE_units",
                   "ARG_MGE_sides_w", "ARG_MGE_net_w", "ARG_MGE_net_w_lo", "ARG_MGE_net_w_hi", "ARG_MGE_net_w_p_enrich", "ARG_MGE_net_w_se_jack",
                   "ARG_MGE_jack_reps"]], on="source", how="left")
    e.columns = [c.replace("ARG_MGE_", "") for c in e.columns]
    t = stats.t.ppf(0.975, (e.jack_reps - 1).clip(lower=1))
    e["net_w_lo_jack"], e["net_w_hi_jack"] = e.net_w - t * e.net_w_se_jack, e.net_w + t * e.net_w_se_jack
    for c in ("their_share", "naive", "net_w"):
        e[f"rank_{c}"] = e[c].rank(ascending=False)
    e = e.sort_values("rank_their_share")
    e.to_csv(f"{out}.env.tsv", sep="\t", index=False)

    rows = []
    ok = e.dropna(subset=["net_w", "naive"])
    for a, b in (("their_share", "naive"), ("their_share", "net_w"), ("naive", "net_w"), ("n50", "naive"), ("n50", "net_w"), ("n50", "their_share")):
        if len(ok) >= 4:
            rho, p = stats.spearmanr(ok[a], ok[b])
            rows.append(dict(level="environment", x=a, y=b, n=len(ok), spearman=rho, p=p))
    m = r[r.fetched == True].merge(tr[["run", "share", "arg_reads"]], left_on="run_accession", right_on="run")  # noqa: E712
    for a, b in (("share", "ARG_MGE_naive"), ("share", "ARG_MGE_net_w"), ("ARG_MGE_naive", "ARG_MGE_net_w"), ("n50", "ARG_MGE_naive"), ("n50", "ARG_MGE_net_w")):
        d = m.dropna(subset=[a, b])
        if len(d) >= 4:
            rho, p = stats.spearmanr(d[a], d[b])
            rows.append(dict(level="run", x=a.replace("ARG_MGE_", "").replace("share", "their_share"), y=b.replace("ARG_MGE_", ""), n=len(d), spearman=rho, p=p))
    rows.append(dict(level="run", x="runs fetched", y="runs with a net_w estimate", n=int((r.fetched == True).sum()),  # noqa: E712
                     spearman=np.nan, p=np.nan, count=int(m.ARG_MGE_net_w.notna().sum())))
    rows.append(dict(level="environment", x="environments", y="net_w interval (package) above 0", n=len(e), spearman=np.nan, p=np.nan,
                     count=int((e.net_w_lo > 0).sum())))
    pd.DataFrame(rows).to_csv(f"{out}.stats.tsv", sep="\t", index=False)
    pd.set_option("display.width", 260)
    pd.set_option("display.max_columns", 40)
    show = ["source", "runs", "projects", "amplified_runs", "n50", "their_share", "naive", "net_w", "net_w_lo", "net_w_hi", "net_w_se_jack",
            "units", "sides_w", "rank_their_share", "rank_naive", "rank_net_w"]
    print(e[show].round(4).to_string(index=False))
    print(pd.DataFrame(rows).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
