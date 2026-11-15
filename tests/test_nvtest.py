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


def test_monte_carlo_null_matches_clump_moments():
    from linkstat.stats.nvtest import mc_counts_mol, mc_counts_w

    rng = np.random.default_rng(3)
    f = rng.random(60) * 0.3
    clump = rng.integers(0, 15, 60)
    sims = mc_counts_w(f, clump, 40_000, np.random.default_rng(4))
    v = sum(clump_var_w(f[clump == k]) for k in np.unique(clump))
    assert abs(sims.mean() - f.sum()) < 0.05 and abs(sims.var() / v - 1) < 0.05
    fl, fr = rng.random(60) * 0.3, rng.random(60) * 0.3
    shared = np.isin(clump, [0, 1, 2])
    sims = mc_counts_mol(fl, fr, clump, shared, 40_000, np.random.default_rng(5))
    p = 1 - (1 - fl) * (1 - fr)
    v = sum(
        clump_var_w(p[clump == k]) if k in (0, 1, 2) else clump_var_mol(fl[clump == k], fr[clump == k])
        for k in np.unique(clump)
    )
    assert abs(sims.mean() - p.sum()) < 0.05 and abs(sims.var() / v - 1) < 0.05


def test_mid_p_values():
    from linkstat.stats.nvtest import mid_p

    assert mid_p(3.0, np.full(99, 3.0)) == (1.0, 0.5, 0.5)
    p, up, low = mid_p(10.0, np.zeros(999))
    assert up == 0.5 / 1000 and p == 2 * up and low > 0.99
    sims = np.random.default_rng(6).poisson(4.0, 100_000).astype(float)
    hits = np.array([mid_p(o, sims)[1] for o in np.random.default_rng(7).poisson(4.0, 4_000)])
    assert abs(np.mean(hits <= 0.05) - 0.05) < 0.015  # mid-p: size near nominal for a discrete count


def test_same_molecule_sides_follow_parity():
    import pandas as pd

    from linkstat.stats.nvtest import score_tests

    curve = (np.array([500.0, 1500.0, 4000.0]), np.array([0.1, 0.2, 0.4]))
    # one locus seen on two reads, the second reversed: its read-left side is the genome-right side
    sd = pd.DataFrame(
        {
            "i": [0, 0, 1, 1],
            "t": [600.0, 3000, 3000, 600],
            "e": [0, 0, 0, 0],
            "side": [0, 1, 0, 1],
            "room": [600.0, 3000, 3000, 600],
        }
    )
    host = np.zeros(4, int)
    unit = sd.i.to_numpy()
    par = score_tests(sd, host, {0: curve}, unit, np.zeros(2, int), np.array([0, 1]), np.zeros(2, bool), 2000, n_mc=0)
    # with parity, both occurrences have the same genome-oriented rooms: they act as one draw (variance 4 q (1 - q))
    q = 1 - (1 - 0.1) * (1 - 0.2)
    assert np.isclose(par["mol"]["V"], 4 * q * (1 - q))
    same = score_tests(sd, host, {0: curve}, unit, np.zeros(2, int), np.zeros(2, int), np.zeros(2, bool), 2000, n_mc=0)
    assert same["mol"]["V"] < par["mol"]["V"]
