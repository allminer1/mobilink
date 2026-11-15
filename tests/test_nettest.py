import numpy as np

from linkstat.simulate import simulate_tables
from linkstat.stats.nettest import net_linked_test

RNG = np.random.default_rng(21)
READS = RNG.lognormal(np.log(9_000), 0.6, 20_000).astype(int) + 500
A_LEN = RNG.integers(600, 1_300, 300)
B_LEN = RNG.integers(600, 1_500, 300)


def _sim(seed, linked, sigma=1.0, n_genomes=100, bases=4e8, a_islands=0.0):
    return simulate_tables(
        READS,
        A_LEN,
        B_LEN,
        seed=seed,
        n_genomes=n_genomes,
        total_bases=bases,
        n_a=24.0 if a_islands == 0 else 12.0,
        a_islands=a_islands,
        islands=5,
        singletons=12,
        linked=linked,
        min_frac=0.8,
        n_truth_anchors=20_000,
        genome_length=(4_000_000, 4_000_000),
        density_sigma=sigma,
        density_cap=2.5,
        positions=True,
    )


def _test(s, strata=True, loci="positions", n_boot=50):
    f, ln = s.features, s.lengths
    return net_linked_test(
        ln,
        f[f["class"] == "ARG"],
        f[f["class"] == "MGE"],
        strata=ln.genome.to_numpy() if strata else None,
        loci=loci,
        n_boot=n_boot,
    )


def test_planted_linkage_is_detected_and_the_interval_covers_truth():
    s = _sim(31, linked=0.3)
    r = _test(s)
    truth = float(s.truth.query("mode == 'end' and w == 2000").net.iloc[0])
    floor = 1 / 10_000  # Monte Carlo reference, n_mc = 9,999 (two-sided mid-p floor)
    assert r["net_w_p"] <= floor and r["net_mol_p"] <= floor and r["net_w_p_enrich"] <= floor / 2
    assert r["net_w_z"] > 5 and r["net_mol_z"] > 5
    assert r["net_w_lo"] - 0.02 < truth < r["net_w_hi"] + 0.02
    assert r["flags"] == "", r["flags"]
    assert r["unit_loci"] < r["units"]  # about 5 reads per locus at this depth


def test_no_linkage_gives_calibrated_scores_and_unit_clumps_matter():
    z, z_noclump = [], []
    for seed in range(40, 52):
        s = _sim(seed, linked=0.0, a_islands=3.0)
        z.append(_test(s)["net_w_z"])
        z_noclump.append(_test(s, loci=None)["net_w_z"])
    assert abs(np.mean(z)) < 0.8 and np.std(z) < 1.6, z
    assert np.std(z_noclump) > 1.5 * np.std(z), (z, z_noclump)  # reads of one locus are not independent


def test_host_strata_remove_host_covariation_under_no_linkage():
    est = {True: [], False: []}
    for seed in range(60, 68):
        s = _sim(seed, linked=0.0, sigma=1.0)
        for strat in (True, False):
            est[strat].append(_test(s, strata=strat, n_boot=20)["net_w"])
    m_strat, m_all = np.mean(est[True]), np.mean(est[False])
    assert abs(m_strat) < 0.004, est[True]
    assert m_all - m_strat > 0.002, (est[True], est[False])


def test_flags_without_strata_or_loci():
    r = _test(_sim(50, linked=0.0, n_genomes=40, bases=1.5e8), strata=False, loci=None, n_boot=20)
    assert "no host strata" in r["flags"] and "no locus clumps" in r["flags"]
    assert np.isfinite(r["net_w_se"])
