"""`mobilink meta`: one estimate per group of runs (e.g. an environment), and contrasts between groups.

The model of the registered survey (amendment 07, Part C):
- The runs of a group are pooled as in `mobilink.pairs.pool_runs`: anchors stay in their run and host stratum.
  The group estimate is the net linked fraction within w of the pooled runs (`net_w_estimate`).
- Runs of one cluster (a project, a site, a study) are not independent. The standard error is the delete-one-cluster
  jackknife: the estimate is recomputed without each cluster in turn, and
  SE = sqrt((k - 1) / k * sum((rep - mean(rep))^2)) over the k replicates with an estimate.
  The interval is estimate +- t(0.975, k - 1) * SE.
- "clusters_effective" is (sum u)^2 / sum(u^2) for the class-A units u of the clusters: the number of equal clusters
  that would carry the same weight. A group that rests on one large cluster has a value near 1.
- Contrasts between two groups: difference of the estimates, SE = sqrt(SE_a^2 + SE_b^2), Welch-Satterthwaite degrees
  of freedom from the jackknife degrees of freedom, t test; Holm adjustment of the two-sided p values over all
  contrasts.
The estimate has no reporting rule here: a replicate must not drop out because it has no event.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from .pairs import pool_runs
from .stats.nettest import net_w_estimate


def jackknife_se(full: float, reps) -> tuple[float, int]:
    """Delete-one-cluster jackknife SE and the number of replicates with an estimate."""
    reps = np.asarray(reps, float)
    reps = reps[np.isfinite(reps)]
    k = len(reps)
    if k < 2 or not np.isfinite(full):
        return np.nan, k
    return float(np.sqrt((k - 1) / k * np.sum((reps - reps.mean()) ** 2))), k


def holm(p) -> np.ndarray:
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.minimum(1, np.maximum.accumulate((len(p) - np.arange(len(p))) * p[o]))
    out = np.empty(len(p))
    out[o] = adj
    return out


def pooled_estimate(runs, a: str, b: str, **kw) -> dict:
    """net_w of class a to class b for a list of run dicts (run, lengths, features, strata); {} if no run or no A call.
    A set of runs in which no A call is left after the same-ORF rule gives a note instead of an estimate."""
    if not runs or not sum(int((r["features"]["class"] == a).sum()) for r in runs):
        return {}
    ln, f, st, _ = pool_runs(runs)
    try:
        return net_w_estimate(ln, f[f["class"] == a], f[f["class"] == b], st, **kw)
    except ValueError as e:
        if "no class A call" not in str(e):
            raise
        return dict(note=str(e))


def group_estimate(runs, clusters, a: str, b: str, **kw) -> tuple[dict, pd.DataFrame]:
    """runs: run dicts; clusters: the cluster label of each run. Returns the group row and the replicate table
    (left_out "none" = all runs)."""
    clusters = [str(c) for c in clusters]
    if len(clusters) != len(runs):
        raise ValueError("one cluster label per run is needed")
    rows = [dict(left_out="none", runs=len(runs), **pooled_estimate(runs, a, b, **kw))]
    for c in sorted(set(clusters)):
        rest = [r for r, x in zip(runs, clusters, strict=True) if x != c]
        rows.append(dict(left_out=c, runs=len(rest), **pooled_estimate(rest, a, b, **kw)))
    rep = pd.DataFrame(rows)
    for c in ("net_w", "units", "sides_w", "events_w"):
        if c not in rep:
            rep[c] = np.nan
    full = rep.iloc[0]
    others = rep.iloc[1:]
    se, k = jackknife_se(full.net_w, others.net_w)
    t = stats.t.ppf(0.975, k - 1) if k >= 2 else np.nan
    u = float(np.nan_to_num(full.units)) - np.nan_to_num(others.units.to_numpy(float))  # A units of each cluster
    eff = float(u.sum() ** 2 / np.sum(u**2)) if np.sum(u**2) > 0 else np.nan
    row = dict(runs=len(runs), clusters=len(set(clusters)), clusters_effective=eff, units=full.units,
               sides_w=full.sides_w, events_w=full.events_w, net_w=full.net_w, se_jack=se, lo95=full.net_w - t * se,
               hi95=full.net_w + t * se, df=k - 1 if k >= 2 else np.nan,
               net_w_min_without_one_cluster=others.net_w.min(), net_w_max_without_one_cluster=others.net_w.max())
    return row, rep


def contrasts(groups: pd.DataFrame, name: str = "group") -> pd.DataFrame:
    """All pairwise contrasts between the groups with a jackknife SE (columns net_w, se_jack, df)."""
    g = groups[np.isfinite(groups.se_jack) & np.isfinite(groups.net_w)].reset_index(drop=True)
    rows = []
    for i in range(len(g)):
        for j in range(i + 1, len(g)):
            ra, rb = g.iloc[i], g.iloc[j]
            d, s = ra.net_w - rb.net_w, float(np.hypot(ra.se_jack, rb.se_jack))
            df = s**4 / (ra.se_jack**4 / ra.df + rb.se_jack**4 / rb.df) if s > 0 else np.nan
            t = d / s if s > 0 else np.nan
            rows.append(dict(a=ra[name], b=rb[name], diff=d, se=s, df=df, t=t, p_a_greater=float(stats.t.sf(t, df)),
                             p_two_sided=float(2 * stats.t.sf(abs(t), df))))
    out = pd.DataFrame(rows, columns=["a", "b", "diff", "se", "df", "t", "p_a_greater", "p_two_sided"])
    if len(out):
        out["p_two_sided_holm"] = holm(out.p_two_sided.to_numpy())
    return out


def meta(manifest: pd.DataFrame, load, a: str, b: str, group: str = "group", cluster: str = "cluster", progress=None,
         **kw) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """manifest: one row per run with the columns run, <group>, <cluster> and whatever `load` needs; load(row) returns
    (lengths, features, strata or None). Groups are processed one after the other, so only the runs of one group are
    in memory. Returns the group table, the replicate table and the contrast table."""
    if manifest.run.duplicated().any():
        raise ValueError("a run is listed twice in the manifest")
    rows, reps = [], []
    for g, sel in manifest.groupby(group, sort=True):
        runs = []
        for r in sel.itertuples(index=False):
            ln, f, st = load(r)
            f = f[f["class"].isin([a, b])][["molecule_id", "start", "end", "class"] + (["scov"] if "scov" in f else [])]
            runs.append(dict(run=r.run, lengths=ln[["molecule_id", "length"]], features=f.reset_index(drop=True),
                             strata=np.zeros(len(ln), int) if st is None else st))
        row, rep = group_estimate(runs, sel[cluster].tolist(), a, b, **kw)
        rows.append({group: g, "class_a": a, "class_b": b, **row})
        reps.append(rep.assign(**{group: g}))
        if progress:
            progress(f"{g}: {row['runs']} runs, {row['clusters']} clusters, net_w {row['net_w']:.4f}")
    groups = pd.DataFrame(rows)
    rep = pd.concat(reps, ignore_index=True)
    rep = rep[[group] + [c for c in rep.columns if c != group]]
    return groups, rep, contrasts(groups, group)
