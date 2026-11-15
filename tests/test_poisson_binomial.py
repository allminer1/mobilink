import itertools

import numpy as np
import pytest

from linkstat.stats import poisson_binomial as pb


def brute_pmf(p):
    out = np.zeros(len(p) + 1)
    for bits in itertools.product((0, 1), repeat=len(p)):
        out[sum(bits)] += np.prod([q if b else 1 - q for q, b in zip(p, bits, strict=True)])
    return out


def test_exact_pmf_matches_enumeration():
    rng = np.random.default_rng(0)
    for n in (1, 2, 5, 11):
        p = rng.random(n)
        assert np.allclose(pb.pmf(p), brute_pmf(p), atol=1e-14)
        for k in range(n + 2):
            assert pb.sf(k, p, method="exact") == pytest.approx(brute_pmf(p)[k:].sum(), abs=1e-14)
            assert pb.cdf(k, p, method="exact") == pytest.approx(brute_pmf(p)[: k + 1].sum(), abs=1e-14)


def test_degenerate_trials_shift_k():
    p = np.array([1.0, 1.0, 0.0, 0.3, 0.6])
    full = brute_pmf(p)
    for k in range(6):
        assert pb.sf(k, p) == pytest.approx(full[k:].sum(), abs=1e-14)
        assert pb.cdf(k, p) == pytest.approx(full[: k + 1].sum(), abs=1e-14)


@pytest.mark.parametrize("n,pmax,seed", [(400, 0.5, 1), (3000, 0.02, 2), (3000, 0.002, 3)])
def test_saddlepoint_tails_have_small_relative_error(n, pmax, seed):
    p = np.random.default_rng(seed).random(n) * pmax
    f = pb.pmf(p)
    mu, var = pb.mean_var(p)
    sd = np.sqrt(var)
    for z in (-4, -2, -1, 0, 1, 2, 4, 6, 8):
        k = int(round(mu + z * sd))
        if k < 0 or k > n:
            continue
        exact_up = f[k:].sum()
        if exact_up > 1e-14 and k >= 1:
            assert pb.sf(k, p, method="saddlepoint") == pytest.approx(exact_up, rel=0.02)
        exact_lo = f[: k + 1].sum()
        if exact_lo > 1e-14 and k < n:
            assert pb.cdf(k, p, method="saddlepoint") == pytest.approx(exact_lo, rel=0.02)


def test_monte_carlo_agrees_with_exact():
    rng = np.random.default_rng(4)
    p = rng.random(200) * 0.1
    draws = (rng.random((100_000, p.size)) < p).sum(axis=1)
    for k in (5, 10, 15):
        mc = (draws >= k).mean()
        se = np.sqrt(mc * (1 - mc) / draws.size)
        assert abs(pb.sf(k, p) - mc) < 4 * se + 1e-4


def test_invalid_probabilities_raise():
    with pytest.raises(ValueError):
        pb.sf(1, [0.2, 1.2])
