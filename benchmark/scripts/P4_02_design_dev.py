"""P4 step 4.2: the meta-analysis design, tried on development data (the 64 public runs of P3; no new run is used).

Purpose: feasibility numbers for amendment 07, not results. From the registered per-run tables of P3:
  1. project-level values: here the random-effects (DerSimonian-Laird) summary of the runs of a project (in P4 the
     runs of a project are pooled with pool_runs; the summary is the cheap stand-in for this look);
  2. class means over projects (random effects) and the contrasts against soil;
  3. meta-regression over runs, with class as fixed effect: naive share on log2(N50), and net_w on log2(N50)
     (random-effects weights for net_w; ordinary least squares for the naive share, which has no standard error).
Output: tables/P4/P4_02_design_dev.projects.tsv, .classes.tsv, .slopes.tsv.
usage: python -I P4_02_design_dev.py
"""
import numpy as np
import pandas as pd
from scipy import stats

L = "/mnt/rmrfc/workflow/iskander/linkstat"


def dl(y, se):
    y, se = np.asarray(y, float), np.asarray(se, float)
    ok = np.isfinite(y) & np.isfinite(se) & (se > 0)
    y, v = y[ok], se[ok] ** 2
    k = len(y)
    if k == 0:
        return dict(k=0, est=np.nan, se=np.nan, tau=np.nan)
    if k == 1:
        return dict(k=1, est=y[0], se=np.sqrt(v[0]), tau=np.nan)
    w = 1 / v
    fixed = np.sum(w * y) / w.sum()
    q = float(np.sum(w * (y - fixed) ** 2))
    c = w.sum() - np.sum(w ** 2) / w.sum()
    tau2 = max(0.0, (q - (k - 1)) / c)
    ws = 1 / (v + tau2)
    return dict(k=k, est=float(np.sum(ws * y) / ws.sum()), se=float(np.sqrt(1 / ws.sum())), tau=float(np.sqrt(tau2)))


def wls(X, y, w):
    XtW = X.T * w
    cov = np.linalg.inv(XtW @ X)
    beta = cov @ XtW @ y
    res = y - X @ beta
    return beta, cov, res


def main():
    pr = pd.concat([pd.read_csv(f"{L}/tables/P3/{d}/per_run.tsv", sep="\t").assign(stage=s) for d, s in (("confirm", 1), ("confirm_ext", 2))])
    pr = pr[(pr.class_a == "ARG") & (pr.class_b == "MGE")]
    proj = pd.concat([pd.read_csv(f"{L}/tables/P3/confirm/p3_final.tsv", sep="\t"), pd.read_csv(f"{L}/public/screen/p3x_final.tsv", sep="\t")])
    pr = pr.merge(proj[["run_accession", "project"]], left_on="run", right_on="run_accession", suffixes=("", "_p"))
    pb = pd.concat([pd.read_csv(f"{L}/tables/P3/{d}/partB.tsv", sep="\t") for d in ("confirm", "confirm_ext")])
    pr = pr.merge(pb[["run", "intact_naive", "intact_n50"]], on="run")
    rows = []
    for (env, p), g in pr.groupby(["env", "project"]):
        r = dl(g.net_w, g.net_w_se)
        rows.append(dict(env=env, project=p, runs=len(g), runs_with_estimate=r["k"], net_w=r["est"], se=r["se"],
                         n50_median=g.intact_n50.median(), naive_median=g.intact_naive.median(), units=int(g.units.sum())))
    P = pd.DataFrame(rows)
    P.to_csv(f"{L}/tables/P4/P4_02_design_dev.projects.tsv", sep="\t", index=False)
    rows = []
    for env, g in P.groupby("env"):
        r = dl(g.net_w, g.se)
        rows.append(dict(env=env, projects=len(g), projects_with_estimate=r["k"], runs=int(g.runs.sum()), net_w=r["est"], se=r["se"], tau_between_projects=r["tau"]))
    C = pd.DataFrame(rows).set_index("env")
    for env in C.index:
        if env != "soil" and np.isfinite(C.se[env]) and np.isfinite(C.se["soil"]):
            d, s = C.net_w[env] - C.net_w["soil"], np.hypot(C.se[env], C.se["soil"])
            C.loc[env, "diff_vs_soil"], C.loc[env, "z"], C.loc[env, "p_one_sided"] = d, d / s, stats.norm.sf(d / s)
    C.to_csv(f"{L}/tables/P4/P4_02_design_dev.classes.tsv", sep="\t")
    # slopes on log2 N50 with class as fixed effect
    rows = []
    envs = sorted(pr.env.unique())
    for name, y, se in (("naive share", pr.intact_naive, None), ("net_w", pr.net_w, pr.net_w_se)):
        ok = np.isfinite(y) & (np.isfinite(se) & (se > 0) if se is not None else True)
        e = pr[ok]
        X = np.column_stack([(e.env == k).to_numpy(float) for k in envs] + [np.log2(e.intact_n50 / 8000)])
        yy = y[ok].to_numpy(float)
        if se is None:
            w = np.ones(len(e))
            beta, cov, res = wls(X, yy, w)
            cov = cov * res.var(ddof=X.shape[1])
        else:
            v = se[ok].to_numpy(float) ** 2
            beta, cov, res = wls(X, yy, 1 / v)
            q = float(np.sum(res ** 2 / v))
            tau2 = max(0.0, (q - (len(yy) - X.shape[1])) / np.sum(1 / v))  # crude moment estimate
            beta, cov, res = wls(X, yy, 1 / (v + tau2))
        b, s = beta[-1], np.sqrt(cov[-1, -1])
        rows.append(dict(outcome=name, runs=len(e), slope_per_doubling_of_N50=b, se=s, lo95=b - 1.96 * s, hi95=b + 1.96 * s,
                         p=2 * stats.norm.sf(abs(b / s)), mean_outcome=float(yy.mean())))
    S = pd.DataFrame(rows)
    S.to_csv(f"{L}/tables/P4/P4_02_design_dev.slopes.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    print(P.round(4).to_string(index=False))
    print(C.round(4).to_string())
    print(S.round(4).to_string(index=False))
    print("runs without an estimate, by class:", pr[~np.isfinite(pr.net_w)].groupby("env").size().to_dict(), "of", pr.groupby("env").size().to_dict())
    print("per-run SE quantiles (net_w):", pr.net_w_se.quantile([.1, .5, .9, .99]).round(4).to_dict())


if __name__ == "__main__":
    main()
