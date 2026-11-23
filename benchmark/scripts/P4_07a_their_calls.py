"""P4 step 4.7, Part 1 (exploratory; D114): the read-level table of Go and Shrestha 2026 (Zenodo 10.5281/zenodo.23492757,
filtered_data.csv: one row per ARG-carrying read, their ARG and MGE calls) and read length.

Their statistic at read level: the share of ARG-carrying reads that also carry an MGE call. MGE call here = insertion
sequence (bounding or standalone), Tn3 / Tn7 / Tn554 marker or integron integrase; plasmid alignments are kept apart.
  (a) the share per run and per environment (pooled reads and mean of the runs);
  (b) share against read length: across runs, and within runs (length classes of the ARG read; logistic slope per
      doubling of read length with run as fixed effect);
  (c) the share after every run's ARG reads are cut in silico to pieces with N50 2 kb (Poisson breakpoints, 5 draws,
      a call cut by a breakpoint stays on a piece if >= 80% of it is there; runs whose ARG reads are shorter are not cut).
Outputs OUT_DIR/their_reads.tsv.gz (one row per ARG read), their_runs.tsv, their_env.tsv, their_length_classes.tsv,
their_slope.tsv.

usage: python -I P4_07a_their_calls.py OUT_DIR      (OUT_DIR: tables/P4/reanalysis)
"""
import ast
import sys

import numpy as np
import pandas as pd
from scipy import stats

L = "/mnt/rmrfc/workflow/iskander/linkstat"
X = f"{L}/external/go_shrestha_2026"
MGE = {"is_bounding": "bounding_IS_coords", "is_standalone": "true_standalone_is_coords", "tn3": "tn3_coords", "tn7": "tn7_coords",
       "tn554": "tn554_coords", "integron": "integron_coords"}
SEED, DRAWS, TARGET = 20261114, 5, 2000
BINS = [0, 1000, 2000, 4000, 8000, 16000, np.inf]


def flat(x):
    """A list of (start, end) pairs from a literal that may be empty, missing or nested one level deeper."""
    if not isinstance(x, str) or not x.strip():
        return []
    v, out = ast.literal_eval(x), []
    stack = list(v) if isinstance(v, (list, tuple)) else []
    while stack:
        y = stack.pop()
        if isinstance(y, (list, tuple)) and len(y) == 2 and all(isinstance(z, (int, float)) for z in y):
            out.append((int(y[0]), int(y[1])))
        elif isinstance(y, (list, tuple)):
            stack.extend(y)
    return out


def load():
    use = ["query_name", "Sample", "query_length", "arg_coords", "plasmid_coords"] + list(MGE.values())
    reads, feats = [], []
    for ch in pd.read_csv(f"{X}/filtered_data.csv", usecols=use, dtype=str, chunksize=200000):
        a = ch.arg_coords.map(flat)
        m = {k: ch[c].map(flat) for k, c in MGE.items()}
        r = pd.DataFrame({"read_id": ch.query_name, "run": ch.Sample, "length": ch.query_length.astype(np.int64), "n_arg": a.map(len),
                          "plasmid": ch.plasmid_coords.map(flat).map(len) > 0, **{k: v.map(len) for k, v in m.items()}})
        reads.append(r)
        for cls, col in [("A", a)] + [("B", v) for v in m.values()]:
            e = col[col.map(len) > 0].explode()
            if len(e):
                feats.append(pd.DataFrame({"read_id": ch.query_name.reindex(e.index).to_numpy(), "start": [p[0] for p in e],
                                           "end": [p[1] for p in e], "cls": cls}))
    r = pd.concat(reads, ignore_index=True)
    r["n_mge"] = r[list(MGE)].sum(axis=1)
    return r, pd.concat(feats, ignore_index=True)


def cut_share(r, f):
    """Share of ARG-carrying pieces with an MGE call after the cut to N50 2 kb; mean of the draws, per run."""
    from mobilink.fragment import fragment_partial, n50, random_bounds, rate_for_n50

    rows = []
    for run, g in r.groupby("run"):
        lengths = g.length.to_numpy(np.int64)
        q = f[f.read_id.isin(set(g.read_id))]
        q = q[(q.end > q.start) & (q.start >= 0)]
        rid = pd.Index(g.read_id).get_indexer(q.read_id)
        s, e = q.start.to_numpy(np.int64), np.minimum(q.end.to_numpy(np.int64), lengths[rid])
        ok = e > s
        d = (lengths, rid[ok], s[ok], e[ok], q.cls.to_numpy(str)[ok])
        native = n50(lengths)
        if native <= TARGET:
            rows.append(dict(run=run, n50_arg_reads=native, cut=False, cut_n50=native, cut_arg_pieces=len(g), cut_share=float((g.n_mge > 0).mean())))
            continue
        rate, v = rate_for_n50(lengths, TARGET, SEED), []
        for k in range(DRAWS):
            pl, piece, _, _, cl = fragment_partial(*d, random_bounds(lengths, rate, SEED + k))
            pa, pb = np.unique(piece[cl == "A"]), np.unique(piece[cl == "B"])
            v.append((n50(pl), len(pa), np.isin(pa, pb).mean() if len(pa) else np.nan))
        v = np.array(v, float)
        rows.append(dict(run=run, n50_arg_reads=native, cut=True, cut_n50=int(v[:, 0].mean()), cut_arg_pieces=int(v[:, 1].mean()),
                         cut_share=float(np.nanmean(v[:, 2]))))
    return pd.DataFrame(rows)


def logit_fe(y, x, group):
    """Logistic regression of y on x with one intercept per group (Newton steps); slope, SE, groups used.
    Groups in which y does not vary carry no information on the slope and are left out."""
    keep = pd.Series(y).groupby(group).transform(lambda s: 0 < s.mean() < 1).to_numpy()
    y, x, group = y[keep], x[keep], np.asarray(group)[keep]
    g = pd.factorize(group)[0]
    Z = np.zeros((len(y), g.max() + 2))
    Z[np.arange(len(y)), g] = 1
    Z[:, -1] = x - x.mean()
    b = np.zeros(Z.shape[1])
    for _ in range(50):
        p = 1 / (1 + np.exp(-Z @ b))
        H = Z.T @ (Z * (p * (1 - p))[:, None])
        step = np.linalg.solve(H, Z.T @ (y - p))
        b += step
        if np.abs(step).max() < 1e-9:
            break
    return float(b[-1]), float(np.sqrt(np.linalg.inv(H)[-1, -1])), int(g.max() + 1), int(len(y))


def main():
    out = sys.argv[1]
    meta = pd.read_csv(f"{L}/public/screen/p4_07_runs.tsv", sep="\t", dtype=str)
    meta["mean_len"] = meta.base_count.astype(float) / meta.read_count.astype(float)
    r, f = load()
    assert set(r.run) <= set(meta.run_accession), sorted(set(r.run) - set(meta.run_accession))[:5]
    r = r.merge(meta[["run_accession", "source", "project_listed"]], left_on="run", right_on="run_accession").drop(columns="run_accession")
    r["has_mge"] = r.n_mge > 0
    print("rows", len(r), "| rows without an ARG call", int((r.n_arg == 0).sum()), "| ARG reads with an MGE call", round(float(r.has_mge.mean()), 4),
          "| with a plasmid alignment", round(float(r.plasmid.mean()), 4), flush=True)
    r.to_csv(f"{out}/their_reads.tsv.gz", sep="\t", index=False)

    c = cut_share(r, f)
    runs = r.groupby("run").agg(source=("source", "first"), project=("project_listed", "first"), arg_reads=("read_id", "size"),
                                median_len_arg_reads=("length", "median"), share=("has_mge", "mean"), share_plasmid=("plasmid", "mean"),
                                **{f"share_{k}": (k, lambda s: float((s > 0).mean())) for k in MGE}).reset_index()
    runs = runs.merge(meta[["run_accession", "read_count", "base_count", "mean_len", "amplified_by_record"]], left_on="run", right_on="run_accession",
                      how="right").drop(columns="run").rename(columns={"run_accession": "run"}).merge(c, on="run", how="left")
    runs["arg_reads_per_read"] = runs.arg_reads / runs.read_count.astype(float)
    runs.to_csv(f"{out}/their_runs.tsv", sep="\t", index=False)

    env = r.groupby("source").agg(runs=("run", "nunique"), projects=("project_listed", "nunique"), arg_reads=("read_id", "size"),
                                  median_len_arg_reads=("length", "median"), share_pooled=("has_mge", "mean")).reset_index()
    rm = runs.dropna(subset=["share"]).groupby("source").agg(share_mean_of_runs=("share", "mean"), cut_share_mean_of_runs=("cut_share", "mean"),
                                                             mean_read_length=("mean_len", "mean"), amplified_runs=("amplified_by_record", lambda s: int((s.astype(str) == "True").sum())))
    w = runs.dropna(subset=["cut_share"])
    pooled_cut = (w.cut_share * w.cut_arg_pieces).groupby(w.source).sum() / w.cut_arg_pieces.groupby(w.source).sum()
    env = env.merge(rm, on="source").merge(pooled_cut.rename("cut_share_pooled").reset_index(), on="source")
    env["rank_native"] = env.share_pooled.rank(ascending=False).astype(int)
    env["rank_cut"] = env.cut_share_pooled.rank(ascending=False).astype(int)
    env.sort_values("rank_native").to_csv(f"{out}/their_env.tsv", sep="\t", index=False)

    r["len_class"] = pd.cut(r.length, BINS, right=False)
    lc = r.groupby(["source", "len_class"], observed=True).agg(arg_reads=("read_id", "size"), share=("has_mge", "mean")).reset_index()
    al = r.groupby("len_class", observed=True).agg(arg_reads=("read_id", "size"), share=("has_mge", "mean")).reset_index().assign(source="all")
    pd.concat([lc, al]).to_csv(f"{out}/their_length_classes.tsv", sep="\t", index=False)

    ok = runs.dropna(subset=["share"])
    rows = []
    for name, x in [("mean read length of the run", ok.mean_len), ("median length of the ARG reads", ok.median_len_arg_reads)]:
        rho, p = stats.spearmanr(x, ok.share)
        rows.append(dict(analysis=f"across runs: Spearman of the share with the {name}", n=len(ok), estimate=rho, p=p))
    b, se, k, n = logit_fe(r.has_mge.to_numpy(float), np.log2(r.length.to_numpy(float)), r.run.to_numpy())
    rows.append(dict(analysis="within runs: log odds of an MGE call per doubling of read length (run fixed effects)", n=n, groups=k, estimate=b, se=se,
                     p=float(2 * stats.norm.sf(abs(b / se)))))
    rho, p = stats.spearmanr(env.share_pooled, env.cut_share_pooled)
    rows.append(dict(analysis="environments: Spearman of the native share with the share after the cut", n=len(env), estimate=rho, p=p))
    pd.DataFrame(rows).to_csv(f"{out}/their_slope.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    print(env.sort_values("rank_native").round(4).to_string(index=False))
    print(pd.DataFrame(rows).round(5).to_string(index=False))
    print(al.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
