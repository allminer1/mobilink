#!/usr/bin/env python3
"""Gate G3b (amendment 04 Part 3), simulation part: false discovery rate from the G3a datasets.

Per setting and test (net_w, net_mol): 1,000 random mixes of 15 datasets under no linkage and 5 at linked share
0.03 (drawn without replacement within a mix), Benjamini-Hochberg at q <= 0.05 on the two-sided p values. False
share of a mix = null datasets among those called (0 if none called). Pass, per test: the mean false share over all
mixes is <= 0.05 + 1.645 x its Monte Carlo SE. Seed 90,500,000.
usage: G3b_summary.py RUN_DIR OUT_TSV
"""
import glob
import sys

import numpy as np
import pandas as pd


def bh(p, q=0.05):
    o = np.argsort(p)
    m = p.size
    ok = p[o] <= q * np.arange(1, m + 1) / m
    k = np.flatnonzero(ok).max() + 1 if ok.any() else 0
    called = np.zeros(m, bool)
    called[o[:k]] = True
    return called


def main():
    run, out = sys.argv[1], sys.argv[2]
    t = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(f"{run}/s*.tsv"))], ignore_index=True)
    rng = np.random.default_rng(90_500_000)
    rows = []
    for k, d in t.groupby("setting"):
        null, lnk = d[d.linked == 0], d[d.linked == 0.03]
        for name in ("net_w", "net_mol"):
            pn, pl = null[f"{name}_p"].to_numpy(), lnk[f"{name}_p"].to_numpy()
            fdr, power = [], []
            for _ in range(1000):
                p = np.r_[rng.choice(pn, 15, replace=False), rng.choice(pl, 5, replace=False)]
                c = bh(p)
                fdr.append(c[:15].sum() / c.sum() if c.any() else 0.0)
                power.append(c[15:].mean())
            rows.append(dict(setting=k, test=name, fdr=np.mean(fdr), fdr_se=np.std(fdr) / np.sqrt(len(fdr)),
                             power=np.mean(power)))
    s = pd.DataFrame(rows)
    s.to_csv(out, sep="\t", index=False)
    for name, g in s.groupby("test"):
        m, se = g.fdr.mean(), np.sqrt(np.sum(g.fdr_se**2)) / len(g)
        print(f"G3b {name}: mean false share {m:.4f} (MC SE {se:.4f}) -> {'pass' if m <= 0.05 + 1.645 * se else 'FAIL'};"
              f" max setting {g.fdr.max():.3f}; mean power {g.power.mean():.2f}")


if __name__ == "__main__":
    main()
