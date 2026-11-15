import numpy as np

from linkstat.simulate import simulate_tables
from linkstat.stats.nettest import net_linked_test

RNG = np.random.default_rng(21)
READS = RNG.lognormal(np.log(9_000), 0.6, 20_000).astype(int) + 500
A_LEN = RNG.integers(600, 1_300, 300)
B_LEN = RNG.integers(600, 1_500, 300)


def _sim(seed, linked, sigma=1.0, n_genomes=100, bases=4e8):
    return simulate_tables(
        READS,
        A_LEN,
        B_LEN,
        seed=seed,
        n_genomes=n_genomes,
        total_bases=bases,
        n_a=24.0,
        islands=5,
        singletons=12,
        linked=linked,
        min_frac=0.8,
        n_truth_anchors=20_000,
        genome_length=(4_000_000, 4_000_000),
        density_sigma=sigma,
    )


def _test(s, strata=True, n_boot=100):
    f, ln = s.features, s.lengths
    return net_linked_test(
        ln,
        f[f["class"] == "ARG"],
        f[f["class"] == "MGE"],
        strata=ln.genome.to_numpy() if strata else None,
        n_boot=n_boot,
    )


def test_planted_linkage_is_detected_and_the_interval_covers_truth():
    s = _sim(31, linked=0.3)
    r = _test(s)
    truth = float(s.truth.query("mode == 'end' and w == 2000").net.iloc[0])
    assert r["net_w_p"] < 1e-6 and r["net_mol_p"] < 1e-6
    assert r["net_w_lo"] - 0.02 < truth < r["net_w_hi"] + 0.02
    assert r["unit"] == "strata" and r["flags"] == ""


def test_host_strata_remove_host_covariation_under_no_linkage():
    est = {True: [], False: []}
    for seed in range(40, 48):
        s = _sim(seed, linked=0.0, sigma=1.0)
        for strat in (True, False):
            est[strat].append(_test(s, strata=strat, n_boot=20)["net_w"])
    m_strat, m_all = np.mean(est[True]), np.mean(est[False])
    assert abs(m_strat) < 0.004, est[True]
    assert m_all - m_strat > 0.002, (est[True], est[False])


def test_without_strata_the_flag_is_set():
    r = _test(_sim(50, linked=0.0, n_genomes=40, bases=1.5e8), strata=False, n_boot=20)
    assert "no host strata" in r["flags"] and r["unit"] == "molecules"


def test_sign_flip_is_calibrated_and_conservative_with_a_dominant_unit():
    from linkstat.stats.nettest import sign_flip

    rng = np.random.default_rng(5)
    p = [sign_flip(rng.normal(0, 1, 40), 999, rng)[0] for _ in range(300)]
    assert 0.02 < np.mean(np.array(p) <= 0.05) < 0.09
    pv = rng.normal(0, 1, 40)
    p0, lo, hi = sign_flip(pv, 1999, rng)
    assert lo < pv.mean() < hi
    assert sign_flip(pv + 3.0, 1999, rng)[0] < 0.001
    # one dominant unit: its sign carries no evidence, so the size stays nominal (a t on G_eff - 1 over-rejects)
    p = [sign_flip(np.r_[50.0 * rng.choice([-1, 1]), rng.normal(0, 0.1, 30)], 499, rng)[0] for _ in range(300)]
    assert np.mean(np.array(p) <= 0.05) < 0.09


def test_expected_events_match_observed_under_no_linkage():
    from linkstat.stats.nettest import expected_events

    rng = np.random.default_rng(6)

    def side(n):
        return dict(t=rng.exponential(20_000, n), e=(rng.random(n) < 0.5).astype(int), w=np.ones(n))

    o, e = expected_events(side(4000), side(20_000), 2000)
    assert abs(o - e) < 4 * np.sqrt(e), (o, e)
