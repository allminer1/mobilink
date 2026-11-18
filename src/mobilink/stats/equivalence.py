"""Equivalence ("at chance") tests on O/E: two one-sided tests (TOST) with bounds 1/r and r (registered: r = 1.25).

H0: O/E >= r or O/E <= 1/r, against H1: 1/r < O/E < r. Each boundary is an alternative of the conditional null with
a common association that shifts the mean of O to r E (or E / r):
- analytic null: the conditional law of O given the margins, tilted (mobilink.stats.conditional, log_tilt) so that its
  saddlepoint mean is the boundary value;
- binned null: in every cell a Fisher noncentral hypergeometric law with one common odds ratio, chosen so that the
  mean of O is the boundary value; the sum over cells is exact (convolution).
p_upper = P(O <= observed | O/E = r), p_lower = P(O >= observed | O/E = 1/r), p_tost = max(p_upper, p_lower).
A boundary outside the support of O (given the margins) cannot hold, so its p value is 0.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import optimize, stats

from . import conditional as cond


def tost_analytic(a, b, n_a, n_b, observed: int, expected: float, bound: float = 1.25, strata=None) -> dict:
    out = {}
    for side, target in (("upper", expected * bound), ("lower", expected / bound)):
        try:
            s = cond.tilt_for_mean(a, b, n_a, n_b, target, strata)
        except ValueError:
            out[f"p_{side}"] = 0.0
            continue
        if side == "upper":
            out["p_upper"] = cond.tail(a, b, n_a, n_b, observed, upper=False, strata=strata, log_tilt=s)
        else:
            out["p_lower"] = cond.tail(a, b, n_a, n_b, observed, upper=True, strata=strata, log_tilt=s)
    out["p_tost"] = max(out["p_upper"], out["p_lower"])
    out["bound"] = bound
    return out


def _nc_pmf(table: pd.DataFrame, psi: float) -> np.ndarray:
    f = np.ones(1)
    for n, a, b in table[["n", "a", "b"]].itertuples(index=False):
        if a == 0 or b == 0:
            continue
        lo, hi = max(0, a + b - n), min(a, b)
        p = stats.nchypergeom_fisher.pmf(np.arange(lo, hi + 1), n, a, b, psi)
        f = np.convolve(f, np.concatenate([np.zeros(lo), p]))
    return f / f.sum()


def _nc_mean(table: pd.DataFrame, psi: float) -> float:
    m = 0.0
    for n, a, b in table[["n", "a", "b"]].itertuples(index=False):
        if a and b:
            m += stats.nchypergeom_fisher.mean(n, a, b, psi)
    return float(m)


def tost_binned(table: pd.DataFrame, observed: int, expected: float, bound: float = 1.25) -> dict:
    """table: cells with columns n, a, b (as mobilink.nulls.binned.cell_table)."""
    t = table[(table.a > 0) & (table.b > 0)]
    lo_s = int(np.sum(np.maximum(0, t.a + t.b - t.n)))
    hi_s = int(np.sum(np.minimum(t.a, t.b)))
    out = {}
    for side, target in (("upper", expected * bound), ("lower", expected / bound)):
        if not lo_s < target < hi_s:
            out[f"p_{side}"] = 0.0
            continue
        lpsi = optimize.brentq(lambda x, tg=target: _nc_mean(t, np.exp(x)) - tg, -30, 30, xtol=1e-10)
        f = _nc_pmf(t, np.exp(lpsi))
        if side == "upper":
            out["p_upper"] = float(f[: observed + 1].sum())
        else:
            out["p_lower"] = float(f[observed:].sum()) if observed < f.size else 0.0
    out["p_tost"] = max(out["p_upper"], out["p_lower"])
    out["bound"] = bound
    return out
