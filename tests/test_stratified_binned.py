import numpy as np
import pandas as pd
import pytest

from mobilink.nulls.analytic import cooccurrence_test
from mobilink.nulls.binned import (
    binned_test,
    cell_table,
    cells,
    length_bins,
    mantel_haenszel,
    null_pmf,
    permutation_test,
)
from mobilink.stats import conditional as cond
from mobilink.stats.equivalence import tost_analytic, tost_binned
from mobilink.stats.multiple import bh


def molecules(rng, n=6_000, lam_a=1 / 40_000, lam_b=1 / 10_000):
    L = np.exp(rng.normal(np.log(6_000), 0.8, n)).astype(int) + 500
    x = np.clip(L - 1_000, 0, None)
    return L, rng.random(n) < -np.expm1(-lam_a * x), rng.random(n) < -np.expm1(-lam_b * x)


def test_exact_binned_null_matches_permutation():
    rng = np.random.default_rng(10)
    L, ha, hb = molecules(rng)
    t = cell_table(cells(length_bins(L)), ha, hb)
    f = null_pmf(t)
    obs, perm = permutation_test(L, ha, hb, n_perm=4_000, seed=1)
    k = np.arange(f.size)
    assert perm.mean() == pytest.approx((f * k).sum(), abs=4 * perm.std() / np.sqrt(perm.size))
    for q in (int(np.quantile(perm, 0.05)), int(np.median(perm)), int(np.quantile(perm, 0.95))):
        mc = (perm >= q).mean()
        assert abs(f[q:].sum() - mc) < 4 * np.sqrt(mc * (1 - mc) / perm.size) + 1e-3
    r = binned_test(L, ha, hb)
    assert r.observed == obs and r.expected == pytest.approx((t.a * t.b / t.n).sum())


def test_stratified_saddlepoint_matches_exact_sum_of_hypergeometrics():
    rng = np.random.default_rng(11)
    sizes = [400, 900, 250, 1200]
    pa, pb_ = [0.05, 0.12, 0.2, 0.08], [0.3, 0.15, 0.4, 0.25]
    g = np.concatenate([np.full(n, s) for s, n in enumerate(sizes)])
    a = np.concatenate([np.full(n, p) for n, p in zip(sizes, pa, strict=True)])
    b = np.concatenate([np.full(n, p) for n, p in zip(sizes, pb_, strict=True)])
    ha, hb = rng.random(g.size) < a, rng.random(g.size) < b
    n_a, n_b = np.bincount(g, ha), np.bincount(g, hb)
    t = pd.DataFrame({"n": sizes, "a": n_a.astype(int), "b": n_b.astype(int)})
    f = null_pmf(t)
    mean = (np.arange(f.size) * f).sum()
    sd = np.sqrt((np.arange(f.size) ** 2 * f).sum() - mean**2)
    for z in (-3, -1.5, 1.5, 3, 4.5):
        k = int(round(mean + z * sd))
        up, lo = f[k:].sum(), f[: k + 1].sum()
        assert cond.tail(a, b, n_a.astype(int), n_b.astype(int), k, True, strata=g) == pytest.approx(up, rel=0.03)
        assert cond.tail(a, b, n_a.astype(int), n_b.astype(int), k, False, strata=g) == pytest.approx(lo, rel=0.03)
    assert cond.variance(a, b, n_a.astype(int), n_b.astype(int), strata=g) == pytest.approx(sd**2, rel=0.03)


def test_mantel_haenszel_single_table_is_woolf():
    t = pd.DataFrame({"n": [100], "a": [30], "b": [40], "ab": [20]})  # 2x2: 20, 10 / 20, 50
    orr, lo, hi = mantel_haenszel(t)
    assert orr == pytest.approx(20 * 50 / (10 * 20))
    se = np.sqrt(1 / 20 + 1 / 10 + 1 / 20 + 1 / 50)
    assert lo == pytest.approx(orr * np.exp(-1.959964 * se), rel=1e-5)
    assert hi == pytest.approx(orr * np.exp(1.959964 * se), rel=1e-5)


def test_bh_known_values_and_nan():
    q = bh([0.01, 0.04, 0.03, 0.005, np.nan])
    assert np.allclose(q[:4], [0.02, 0.04, 0.04, 0.02]) and np.isnan(q[4])


def test_tost_analytic_equals_binned_for_one_cell_with_equal_lengths():
    rng = np.random.default_rng(12)
    n = 4_000
    ha, hb = rng.random(n) < 0.1, rng.random(n) < 0.3
    a, b = np.full(n, ha.mean()), np.full(n, hb.mean())
    obs, e = int((ha & hb).sum()), float(ha.sum() * hb.sum() / n)
    t = pd.DataFrame({"n": [n], "a": [int(ha.sum())], "b": [int(hb.sum())]})
    ta = tost_analytic(a, b, int(ha.sum()), int(hb.sum()), obs, e)
    tb = tost_binned(t, obs, e)
    assert ta["p_upper"] == pytest.approx(tb["p_upper"], rel=0.05, abs=1e-6)
    assert ta["p_lower"] == pytest.approx(tb["p_lower"], rel=0.05, abs=1e-6)
    assert ta["p_tost"] < 0.05  # 4,000 independent molecules: equivalent to chance within 1.25


def test_tost_is_not_significant_under_association():
    rng = np.random.default_rng(13)
    L, ha, hb = molecules(rng, lam_b=1 / 60_000)  # B rare, so a planted association moves O/E well past 1.25
    hb = hb | (ha & (rng.random(L.size) < 0.3))
    r = cooccurrence_test(L, ha, hb, np.array([1_000.0]), np.array([1_000.0]), min_frac=1.0)
    assert r.oe > 1.5 and r.p_tost > 0.5


def test_stratified_null_removes_confounding_by_sample():
    """Two samples with different densities of both classes: pooled, the length-only null sees spurious enrichment;
    fitted per sample (N1), O/E returns to about 1."""
    rng = np.random.default_rng(14)
    parts = [molecules(rng, 8_000, 1 / 80_000, 1 / 30_000), molecules(rng, 8_000, 1 / 8_000, 1 / 3_000)]
    L = np.concatenate([p[0] for p in parts])
    ha = np.concatenate([p[1] for p in parts])
    hb = np.concatenate([p[2] for p in parts])
    strata = np.repeat(["s1", "s2"], 8_000)
    pooled = cooccurrence_test(L, ha, hb, np.array([1_000.0]), np.array([1_000.0]), min_frac=1.0, tost_bound=None)
    strat = cooccurrence_test(
        L, ha, hb, np.array([1_000.0]), np.array([1_000.0]), strata=strata, min_frac=1.0, tost_bound=None
    )
    assert pooled.oe > 1.15 and pooled.p_enriched < 1e-4
    assert abs(strat.oe - 1) < 0.08 and strat.p_two_sided > 0.01 and strat.n_strata == 2
    b = binned_test(L, ha, hb, strata=strata)
    assert abs(b.oe - 1) < 0.1
