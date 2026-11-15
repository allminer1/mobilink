import numpy as np

from linkstat.stats.nvtest import cdf_at, clump_var_mol, clump_var_w, km_curve


def test_clump_variances_match_their_definitions():
    f = 0.3
    assert np.isclose(clump_var_w(np.array([f])), f * (1 - f))
    assert np.isclose(clump_var_w(np.full(4, f)), 16 * f * (1 - f))  # one shared draw seen 4 times
    rng = np.random.default_rng(1)
    fv = rng.random(5)
    d = rng.random(200_000)  # D ~ U(0, 1), so F(tau) = tau
    y = (d[:, None] <= fv[None, :]).sum(1)
    assert abs(y.var() - clump_var_w(fv)) < 0.03
    fl, fr = rng.random(4) * 0.5, rng.random(4) * 0.5
    dl, dr = rng.random(200_000), rng.random(200_000)
    ym = ((dl[:, None] <= fl[None, :]) | (dr[:, None] <= fr[None, :])).sum(1)
    assert abs(ym.var() - clump_var_mol(fl, fr)) < 0.03


def test_km_curve_and_cdf():
    t = np.array([1.0, 2.0, 3.0, 4.0])
    e = np.array([1, 0, 1, 0])
    c = km_curve(t, e, np.ones(4))
    assert np.allclose(cdf_at(c, [0.5, 1, 2.5, 3, 10]), [0, 0.25, 0.25, 1 - 0.75 * 0.5, 1 - 0.75 * 0.5])
    assert np.allclose(cdf_at((np.array([]), np.array([])), [1, 2]), 0)
