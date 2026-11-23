"""P4 step 4.7b (exploratory; D115): one sample sequenced after several extraction or library methods.

Input: the cell table of P4_07c_estimates.py for public/screen/p4_07b_runs.tsv, grouped by set,sample,condition
(cell = sample x method, replicate runs pooled).
  slopes     weighted least squares of the naive share and of the net linked fraction (ARG:MGE; MRG:MGE as a second
             pair) on log2 read N50 of the cell, one intercept per sample, weights = A-carrying reads (naive share) or
             A units (net linked fraction) of the cell; cluster-robust standard error by sample; per set and over all
             sets. Cells without an estimate are left out of that outcome; samples with one cell left are dropped.
  profiles   Bray-Curtis dissimilarity of the Kraken2 genus profiles (reads at confidence 0, runs of a cell summed):
             between methods of one sample, and between samples of one set under the same method. It shows how far
             "the same community" holds between methods.
Outputs OUT_PREFIX.slopes.tsv, OUT_PREFIX.cells.tsv (cells with the centred values used), OUT_PREFIX.profiles.tsv.

usage: python -I P4_07d_kit_slopes.py CELLS_TSV RUNS_TSV OUT_PREFIX
"""
import sys

import numpy as np
import pandas as pd
from scipy import stats

L = "/mnt/rmrfc/workflow/iskander/linkstat"


def wslope(d, y, w):
    d = d[np.isfinite(d[y]) & np.isfinite(d.x) & (d[w] > 0)].copy()
    d = d[d.groupby("g")[y].transform("size") >= 2]
    n, G = len(d), d.g.nunique()
    if n < 4 or G < 2:
        return dict(cells=n, samples=G), d
    ww = d[w] / d[w].sum()
    mx, my = (ww * d.x).groupby(d.g).transform("sum") / ww.groupby(d.g).transform("sum"), (ww * d[y]).groupby(d.g).transform("sum") / ww.groupby(d.g).transform("sum")
    xc, yc = d.x - mx, d[y] - my
    sxx = float((ww * xc ** 2).sum())
    b = float((ww * xc * yc).sum() / sxx)
    e = yc - b * xc
    meat = float((((ww * xc * e).groupby(d.g).sum()) ** 2).sum())
    se = float(np.sqrt(meat / sxx ** 2 * G / (G - 1) * (n - 1) / max(n - G - 1, 1)))
    t90 = stats.t.ppf(0.95, G - 1)
    d[f"{y}_xc"], d[f"{y}_yc"] = xc, yc
    return dict(cells=n, samples=G, slope=b, se=se, lo90=b - t90 * se, hi90=b + t90 * se, p_two_sided=float(2 * stats.t.sf(abs(b / se), G - 1)),
                n50_min=float(2 ** d.x.min()), n50_max=float(2 ** d.x.max()),
                largest_n50_ratio_within_sample=float(2 ** (d.x.groupby(d.g).max() - d.x.groupby(d.g).min()).max())), d


def genus_profile(runs):
    tot = None
    for r in runs:
        t = pd.read_csv(f"{L}/public/{r}/{r}.kraken_c0.report", sep="\t", header=None, usecols=[1, 3, 5], names=["reads", "rank", "name"])
        g = t[t["rank"] == "G"].assign(name=lambda d: d.name.str.strip()).groupby("name").reads.sum()
        tot = g if tot is None else tot.add(g, fill_value=0)
    return tot / tot.sum()


def bray(a, b):
    i = a.index.union(b.index)
    a, b = a.reindex(i, fill_value=0), b.reindex(i, fill_value=0)
    return float(1 - np.minimum(a, b).sum())


def main():
    cells = pd.read_csv(sys.argv[1], sep="\t")
    runs = pd.read_csv(sys.argv[2], sep="\t", dtype=str)
    out = sys.argv[3]
    cells = cells[cells.runs > 0].assign(x=lambda d: np.log2(d.n50), g=lambda d: d["set"] + "|" + d["sample"])
    rows, keep = [], cells.copy()
    for pair in ("ARG_MGE", "MRG_MGE"):
        for y, w, name in ((f"{pair}_naive", f"{pair}_reads_A", "naive share"), (f"{pair}_net_w", f"{pair}_units", "net linked fraction within 2 kb")):
            if y not in cells:
                continue
            for s in list(cells["set"].unique()) + ["all sets"]:
                d = cells if s == "all sets" else cells[cells["set"] == s]
                r, used = wslope(d, y, w)
                rows.append(dict(pair=pair.replace("_", ":"), outcome=name, set=s, **r))
                if s == "all sets" and f"{y}_xc" in used:
                    keep = keep.merge(used[["set", "sample", "condition", f"{y}_xc", f"{y}_yc"]], on=["set", "sample", "condition"], how="left")
    pd.DataFrame(rows).to_csv(f"{out}.slopes.tsv", sep="\t", index=False)
    keep.to_csv(f"{out}.cells.tsv", sep="\t", index=False)

    import os
    have = runs[[os.path.exists(f"{L}/public/{r}/{r}.kraken_c0.report") for r in runs.run_accession]]
    prof = {k: genus_profile(g.run_accession) for k, g in have.groupby(["set", "sample", "condition"])}
    pr = []
    keys = list(prof)
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            if a[0] != b[0]:
                continue
            kind = "same sample, other method" if a[1] == b[1] else "other sample, same method" if a[2] == b[2] else None
            if kind:
                pr.append(dict(set=a[0], kind=kind, a="|".join(a[1:]), b="|".join(b[1:]), bray_curtis=bray(prof[a], prof[b])))
    pr = pd.DataFrame(pr)
    pr.to_csv(f"{out}.profiles.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    print(pd.DataFrame(rows).round(4).to_string(index=False))
    if len(pr):
        print(pr.groupby(["set", "kind"]).bray_curtis.describe()[["count", "50%", "min", "max"]].round(3).to_string())


if __name__ == "__main__":
    main()
