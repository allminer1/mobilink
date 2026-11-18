"""N2: the permutation null within length bins (and optional strata), computed exactly.

Shuffling the B labels among the molecules of each cell (length bin x stratum) keeps the cell margins. The number of
molecules with A and B in a cell is then hypergeometric (n_c molecules, a_c with A, b_c with B), and the total O is a
sum of independent hypergeometric counts. Its exact distribution is computed by convolution, so no Monte Carlo is
needed; this is also the exact conditional (CMH-type) test across the cells. The expectation is
E = sum_c a_c b_c / n_c, the prototype's "binned" null. `permutation_test` draws the same null by Monte Carlo, as a
check and for statistics without a closed form.
Length bins follow the prototype: width 1 kb, molecules above the 99.5% length quantile pooled into the last bin.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd
from scipy import stats

from .analytic import poisson_interval


def length_bins(lengths, width: int = 1000, cap_quantile: float = 0.995) -> np.ndarray:
    L = np.asarray(lengths, dtype=float)
    cap = np.quantile(L, cap_quantile)
    return (np.minimum(L, cap) // width).astype(np.int64)


def cells(bins, strata=None) -> np.ndarray:
    """Integer cell code per molecule for length bin x stratum."""
    key = pd.Series(bins) if strata is None else pd.Series(list(zip(strata, bins, strict=True)))
    return pd.factorize(key)[0].astype(np.int64)


def cell_table(cell, has_a, has_b) -> pd.DataFrame:
    d = pd.DataFrame({"c": cell, "a": np.asarray(has_a, bool), "b": np.asarray(has_b, bool)})
    d["ab"] = d.a & d.b
    return d.groupby("c").agg(n=("a", "size"), a=("a", "sum"), b=("b", "sum"), ab=("ab", "sum"))


def null_pmf(table: pd.DataFrame) -> np.ndarray:
    """Exact pmf of O = sum of independent hypergeometric counts over the cells."""
    f = np.ones(1)
    for n, a, b in table[["n", "a", "b"]].itertuples(index=False):
        if a == 0 or b == 0:
            continue
        lo, hi = max(0, a + b - n), min(a, b)
        p = stats.hypergeom.pmf(np.arange(lo, hi + 1), n, a, b)
        f = np.convolve(f, np.concatenate([np.zeros(lo), p]))
    return f / f.sum()


@dataclass
class BinnedResult:
    n_molecules: int
    n_cells: int
    observed: int
    expected: float
    var_null: float
    oe: float
    oe_lo: float
    oe_hi: float
    p_enriched: float
    p_depleted: float
    p_two_sided: float
    p_enriched_mid: float
    p_depleted_mid: float
    p_two_sided_mid: float
    mh_odds_ratio: float  # Mantel-Haenszel common odds ratio over the cells
    mh_or_lo: float  # 95% interval (Robins-Breslow-Greenland)
    mh_or_hi: float

    def as_dict(self):
        return asdict(self)


def mantel_haenszel(table: pd.DataFrame, alpha: float = 0.05):
    """Common odds ratio over the cells with the Robins-Breslow-Greenland variance of its log."""
    n = table.n.to_numpy(float)
    a11 = table.ab.to_numpy(float)
    a10 = table.a.to_numpy(float) - a11
    a01 = table.b.to_numpy(float) - a11
    a00 = n - a11 - a10 - a01
    ok = n > 1
    n, a11, a10, a01, a00 = n[ok], a11[ok], a10[ok], a01[ok], a00[ok]
    r, s = a11 * a00 / n, a10 * a01 / n
    R, S = r.sum(), s.sum()
    if R == 0 or S == 0:
        return (np.inf if S == 0 and R > 0 else 0.0 if R == 0 and S > 0 else np.nan), np.nan, np.nan
    orr = R / S
    p, q = (a11 + a00) / n, (a10 + a01) / n
    var = (p * r).sum() / (2 * R**2) + ((p * s).sum() + (q * r).sum()) / (2 * R * S) + (q * s).sum() / (2 * S**2)
    zq = stats.norm.ppf(1 - alpha / 2)
    return float(orr), float(orr * np.exp(-zq * np.sqrt(var))), float(orr * np.exp(zq * np.sqrt(var)))


def binned_test(
    lengths,
    has_a,
    has_b,
    strata=None,
    width: int = 1000,
    cap_quantile: float = 0.995,
    alpha: float = 0.05,
    bin_by=None,
) -> BinnedResult:
    # bin_by: the length used for the bins (default the molecule length). With bin_by = length - bases covered by B
    # (mobilink.features.covered_bases), molecules in a cell have the same room for A, which removes the
    # excluded-volume bias of length-only bins.
    cell = cells(length_bins(lengths if bin_by is None else bin_by, width, cap_quantile), strata)
    t = cell_table(cell, has_a, has_b)
    obs = int(t.ab.sum())
    n, a, b = t.n.to_numpy(float), t.a.to_numpy(float), t.b.to_numpy(float)
    mu = float(np.sum(a * b / n))
    var = float(np.sum(np.where(n > 1, a * b * (n - a) * (n - b) / (n**2 * np.maximum(n - 1, 1)), 0.0)))
    f = null_pmf(t)  # f[k] = P(O = k)
    k = obs
    up = float(f[k:].sum()) if k < f.size else 0.0
    dn = float(f[: k + 1].sum())
    pk = float(f[k]) if k < f.size else 0.0
    lo, hi = poisson_interval(obs, alpha)
    orr, or_lo, or_hi = mantel_haenszel(t, alpha)
    return BinnedResult(
        n_molecules=int(n.sum()),
        n_cells=int(len(t)),
        observed=obs,
        expected=mu,
        var_null=var,
        oe=obs / mu if mu > 0 else np.nan,
        oe_lo=lo / mu if mu > 0 else np.nan,
        oe_hi=hi / mu if mu > 0 else np.nan,
        p_enriched=min(1.0, up),
        p_depleted=min(1.0, dn),
        p_two_sided=min(1.0, 2 * min(up, dn)),
        p_enriched_mid=up - pk / 2,
        p_depleted_mid=dn - pk / 2,
        p_two_sided_mid=min(1.0, 2 * min(up, dn) - pk),
        mh_odds_ratio=orr,
        mh_or_lo=or_lo,
        mh_or_hi=or_hi,
    )


def permutation_test(
    lengths,
    has_a,
    has_b,
    strata=None,
    n_perm: int = 2000,
    seed: int = 20261114,
    width: int = 1000,
    cap_quantile: float = 0.995,
):
    """Monte Carlo version of the same null: B labels shuffled within cells. Returns (observed, null counts)."""
    rng = np.random.default_rng(seed)
    cell = cells(length_bins(lengths, width, cap_quantile), strata)
    has_a, has_b = np.asarray(has_a, bool), np.asarray(has_b, bool)
    order = np.argsort(cell, kind="stable")
    bounds = np.flatnonzero(np.diff(cell[order])) + 1
    groups = np.split(order, bounds)
    out = np.empty(n_perm, dtype=np.int64)
    a_sorted = [has_a[g] for g in groups]
    b_sorted = [has_b[g] for g in groups]
    for i in range(n_perm):
        out[i] = sum(int(np.sum(a & rng.permutation(b))) for a, b in zip(a_sorted, b_sorted, strict=True))
    return int(np.sum(has_a & has_b)), out
