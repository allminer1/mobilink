import numpy as np
import pandas as pd
import pytest

from linkstat.features import covered_bases
from linkstat.nulls.analytic import clustering_index, cooccurrence_test, excluded_volume_exposure, exposure
from linkstat.nulls.binned import binned_test
from linkstat.simulate import simulate_tables


def test_covered_bases_merges_overlaps():
    f = pd.DataFrame(
        {
            "molecule_id": ["m1", "m1", "m1", "m2", "m3", "m3"],
            "start": [0, 50, 300, 10, 0, 100],
            "end": [100, 120, 400, 20, 500, 200],
        }
    )
    assert list(covered_bases(f, ["m1", "m2", "m3", "m4"])) == [220.0, 10.0, 500.0, 0.0]


def test_excluded_volume_exposure_by_hand():
    # molecule r1 (10 kb): B at [4000, 5000) and [5200, 6000); A calls of 900 bp; f = 0.8
    #   head gap 4000 -> 4000 - 0.8 * 900 = 3280; interior gap 200 -> 0; tail gap 4000 -> 3280
    # molecule r2 (10 kb) without B: 10000 - 0.6 * 900 = 9460
    b = pd.DataFrame({"molecule_id": ["r1", "r1"], "start": [4000, 5200], "end": [5000, 6000]})
    x = excluded_volume_exposure([10_000, 10_000], ["r1", "r2"], b, np.array([900.0]), 0.8)
    assert x == pytest.approx([6560.0, 9460.0])
    assert np.allclose(
        excluded_volume_exposure([5_000, 9_000], ["r1", "r2"], b.iloc[:0], np.array([900.0])),
        exposure([5_000, 9_000], np.array([900.0])),
    )


def test_excluded_volume_correction_removes_the_depletion_bias():
    rng = np.random.default_rng(0)
    rl = rng.lognormal(np.log(7_000), 0.6, 20_000).astype(int) + 500
    al = np.maximum(150, rng.lognormal(np.log(900), 0.35, 2_000)).astype(int)
    bl = np.maximum(150, rng.lognormal(np.log(1000), 0.5, 2_000)).astype(int)
    oe = {"N0": [], "N0ev": [], "N2": [], "N2ev": []}
    for seed in range(4):
        sim = simulate_tables(
            rl,
            al,
            bl,
            seed=200 + seed,
            n_genomes=200,
            total_bases=1.6e8,
            n_a=24.0,
            islands=5,
            singletons=12,
            linked=0.0,
            min_frac=0.8,
            n_truth_anchors=1_000,
            genome_length=(4_000_000, 4_000_000),
            abundance_sigma=0.0,
        )
        f, ids, ln = sim.features, sim.lengths.molecule_id, sim.lengths.length.to_numpy()
        fa, fb = f[f["class"] == "ARG"], f[f["class"] == "MGE"]
        ha, hb = ids.isin(fa.molecule_id).to_numpy(), ids.isin(fb.molecule_id).to_numpy()
        la, lb = (fa.end - fa.start).to_numpy(float), (fb.end - fb.start).to_numpy(float)
        xa = excluded_volume_exposure(ln, ids, fb, la, 0.8)
        oe["N0"].append(cooccurrence_test(ln, ha, hb, la, lb, tost_bound=None).oe)
        r = cooccurrence_test(ln, ha, hb, la, lb, tost_bound=None, exposure_a=xa)
        assert r.excluded_volume
        oe["N0ev"].append(r.oe)
        oe["N2"].append(binned_test(ln, ha, hb).oe)
        oe["N2ev"].append(binned_test(ln, ha, hb, bin_by=xa).oe)
    m = {k: np.mean(v) for k, v in oe.items()}
    assert m["N0"] < 0.85 and m["N2"] < 0.9  # length-only nulls: biased low by excluded volume
    assert abs(m["N0ev"] - 1) < 0.06 and abs(m["N2ev"] - 1) < 0.08


def test_ev_helper_conditions_on_the_more_common_class():
    from linkstat.nulls.analytic import cooccurrence_test_ev

    rng = np.random.default_rng(5)
    rl = rng.lognormal(np.log(7_000), 0.6, 10_000).astype(int) + 500
    al = np.maximum(150, rng.lognormal(np.log(900), 0.35, 1_000)).astype(int)
    bl = np.maximum(150, rng.lognormal(np.log(1000), 0.5, 1_000)).astype(int)
    sim = simulate_tables(
        rl,
        al,
        bl,
        seed=9,
        n_genomes=100,
        total_bases=1e8,
        n_a=24.0,
        islands=5,
        singletons=12,
        linked=0.0,
        n_truth_anchors=1_000,
    )
    f = sim.features
    fa, fb = f[f["class"] == "ARG"], f[f["class"] == "MGE"]
    r1, c1, _ = cooccurrence_test_ev(sim.lengths, fa, fb, tost_bound=None)
    r2, c2, _ = cooccurrence_test_ev(sim.lengths, fb, fa, tost_bound=None)
    assert r1.observed == r2.observed and r1.expected == pytest.approx(r2.expected)
    assert {c1, c2} == {"A", "B"} and r1.excluded_volume


def test_clustering_index_is_one_for_uniform_genes_and_large_for_islands():
    rng = np.random.default_rng(11)
    reads = rng.lognormal(np.log(12_000), 0.6, 20_000).astype(int) + 500
    s = simulate_tables(
        reads,
        rng.integers(600, 1_400, 300),
        rng.integers(600, 1_400, 300),
        seed=5,
        n_genomes=60,
        total_bases=1e9,
        n_a=24.0,
        islands=5,
        singletons=12,
        clustered=True,
        linked=0.0,
        min_frac=0.8,
        n_truth_anchors=1000,
    )
    f, ln = s.features, s.lengths
    ia = clustering_index(ln.length.to_numpy(), ln.molecule_id, f[f["class"] == "ARG"])
    ib = clustering_index(ln.length.to_numpy(), ln.molecule_id, f[f["class"] == "MGE"])
    assert 0.8 < ia < 1.25 and ib > 3
