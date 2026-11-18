import numpy as np
import pytest

from mobilink.nulls.analytic import cooccurrence_test, exposure, fit_density, poisson_interval, presence


def test_exposure_matches_brute_force_mixture():
    rng = np.random.default_rng(0)
    L = rng.integers(0, 20_000, 300)
    a = rng.integers(300, 3_000, 50)
    w = rng.random(50)
    for f in (1.0, 0.8, 0.5):
        c = 2 * f - 1
        brute = np.array([np.sum(w / w.sum() * np.clip(x - c * a, 0, None)) for x in L])
        assert np.allclose(exposure(L, a, f, w), brute)
    assert np.allclose(exposure(L, 900, 1.0), np.clip(L - 900, 0, None))


def test_fit_density_reproduces_count():
    x = np.random.default_rng(1).integers(0, 30_000, 5_000).astype(float)
    lam = fit_density(x, 700)
    assert presence(x, lam).sum() == pytest.approx(700, rel=1e-10)
    assert fit_density(x, 0) == 0.0


def test_poisson_interval_known_values():
    lo, hi = poisson_interval(0)
    assert lo == 0 and hi == pytest.approx(3.689, abs=1e-3)
    lo, hi = poisson_interval(10)
    assert lo == pytest.approx(4.795, abs=1e-3) and hi == pytest.approx(18.390, abs=1e-3)


def independent_molecules(rng, n=20_000, lam_a=1 / 30_000, lam_b=1 / 15_000, gene=1_000):
    """Molecules with log-normal lengths; A and B present independently given length (exposure model)."""
    L = np.exp(rng.normal(np.log(6_000), 0.8, n)).astype(int)
    x = np.clip(L - gene, 0, None)
    has_a = rng.random(n) < -np.expm1(-lam_a * x)
    has_b = rng.random(n) < -np.expm1(-lam_b * x)
    return L, has_a, has_b


def test_oe_near_one_and_conditional_test_has_nominal_size():
    rng = np.random.default_rng(2)
    p_cond, p_unc, oes = [], [], []
    for _ in range(400):
        L, ha, hb = independent_molecules(rng)
        r = cooccurrence_test(L, ha, hb, np.array([1_000.0]), np.array([1_000.0]), min_frac=1.0)
        p_cond.append(r.p_enriched)
        p_unc.append(r.p_enriched_unconditional)
        oes.append(r.oe)
    assert abs(np.mean(oes) - 1) < 0.03
    # one-sided test at 5%; 400 datasets: 95% band about 2.9-7.3%
    assert 0.029 <= np.mean(np.array(p_cond) <= 0.05) <= 0.073
    # the plug-in Poisson-binomial test is conservative when the classes are common (here about 20% and 35%)
    assert np.mean(np.array(p_unc) <= 0.05) < 0.03


def test_conditional_variance_is_hypergeometric_at_equal_lengths():
    rng = np.random.default_rng(6)
    n = 5_000
    L = np.full(n, 8_000)
    ha, hb = rng.random(n) < 0.2, rng.random(n) < 0.3
    r = cooccurrence_test(L, ha, hb, np.array([1_000.0]), np.array([1_000.0]), min_frac=1.0)
    a, b = ha.mean(), hb.mean()
    assert r.var_null == pytest.approx(n * a * b * (1 - a) * (1 - b), rel=1e-9)
    assert r.var_unconditional == pytest.approx(n * a * b * (1 - a * b), rel=1e-9)


def test_planted_linkage_is_detected():
    rng = np.random.default_rng(3)
    L, ha, hb = independent_molecules(rng)
    hb = hb | (ha & (rng.random(L.size) < 0.2))
    r = cooccurrence_test(L, ha, hb, np.array([1_000.0]), np.array([1_000.0]), min_frac=1.0)
    # refitting lambda_B absorbs part of the planted signal, so O/E is modest, but the test is decisive
    assert r.oe > 1.1 and r.oe_lo > 1 and r.p_enriched < 1e-10


def test_mean_mode_uses_whole_genes():
    rng = np.random.default_rng(4)
    L, ha, hb = independent_molecules(rng, n=5_000)
    r = cooccurrence_test(L, ha, hb, np.array([800.0, 1_200.0]), np.array([1_000.0]), mode="mean")
    assert r.min_frac == 1.0 and r.exposure_mode == "mean"


def test_mid_p_lies_between_the_two_tails():
    rng = np.random.default_rng(9)
    L, ha, hb = independent_molecules(rng, n=3_000, lam_a=1 / 300_000, lam_b=1 / 150_000)
    r = cooccurrence_test(L, ha, hb, np.array([1_000.0]), np.array([1_000.0]), min_frac=1.0)
    assert r.p_enriched_mid <= r.p_enriched and r.p_depleted_mid <= r.p_depleted
    assert r.p_enriched_mid + r.p_depleted_mid == pytest.approx(1.0, abs=0.02)
