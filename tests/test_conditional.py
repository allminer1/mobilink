import numpy as np
import pytest
from scipy import stats

from linkstat.nulls.analytic import conditional_variance
from linkstat.stats import conditional as cond


def test_exact_pmf_is_hypergeometric_for_equal_probabilities():
    n, na, nb = 60, 12, 20
    a, b = np.full(n, 0.2), np.full(n, 0.35)
    f = cond.exact_pmf(a, b, na, nb)
    ref = stats.hypergeom.pmf(np.arange(min(na, nb) + 1), n, na, nb)
    assert np.allclose(f, ref, atol=1e-12)


@pytest.mark.parametrize("n,na,nb", [(2000, 200, 600), (5000, 50, 400), (800, 30, 30)])
def test_saddlepoint_tail_matches_hypergeometric(n, na, nb):
    a, b = np.full(n, na / n), np.full(n, nb / n)
    mean = na * nb / n
    sd = np.sqrt(stats.hypergeom.var(n, na, nb))
    for zz in (-4, -2, -1, 1, 2, 4, 6):
        k = max(0, int(round(mean + zz * sd)))
        if k > min(na, nb):
            continue
        up = stats.hypergeom.sf(k - 1, n, na, nb)
        lo = stats.hypergeom.cdf(k, n, na, nb)
        if up > 1e-12:
            assert cond.tail(a, b, na, nb, k, upper=True) == pytest.approx(up, rel=0.03)
        if lo > 1e-12:
            # P(O <= 0) sits on the edge of the support, where the continuity-corrected approximation is ~5% low
            assert cond.tail(a, b, na, nb, k, upper=False) == pytest.approx(lo, rel=0.03 if k > 0 else 0.06)


def test_saddlepoint_tail_matches_exact_recursion_with_unequal_lengths():
    rng = np.random.default_rng(7)
    n = 250
    L = np.exp(rng.normal(np.log(6000), 0.9, n))
    x = np.clip(L - 1000, 0, None)
    a, b = -np.expm1(-x / 40_000), -np.expm1(-x / 12_000)
    na, nb = int(round(a.sum())), int(round(b.sum()))
    f = cond.exact_pmf(a, b, na, nb)
    for k in range(len(f)):
        up, lo = f[k:].sum(), f[: k + 1].sum()
        if 1e-10 < up < 0.999 and k > 0:
            assert cond.tail(a, b, na, nb, k, upper=True) == pytest.approx(up, rel=0.05)
        if 1e-10 < lo < 0.999:
            assert cond.tail(a, b, na, nb, k, upper=False) == pytest.approx(lo, rel=0.05)
    # the conditional variance at the saddlepoint matches the exact one
    o = np.arange(len(f))
    exact_var = (f * o**2).sum() - (f * o).sum() ** 2
    assert cond.variance(a, b, na, nb) == pytest.approx(exact_var, rel=0.05)


def test_saddlepoint_variance_agrees_with_delta_method():
    rng = np.random.default_rng(8)
    n = 20_000
    L = np.exp(rng.normal(np.log(6000), 0.8, n))
    x = np.clip(L - 1000, 0, None)
    la, lb = 1 / 30_000, 1 / 15_000
    a, b = -np.expm1(-la * x), -np.expm1(-lb * x)
    na, nb = int(round(a.sum())), int(round(b.sum()))
    assert cond.variance(a, b, na, nb) == pytest.approx(conditional_variance(x, x, la, lb), rel=0.03)
