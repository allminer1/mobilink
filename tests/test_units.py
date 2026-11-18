import numpy as np
import pandas as pd

from mobilink.nulls.analytic import excluded_volume_exposure
from mobilink.nulls.units import Catalogue, cooccurrence_test_units, make_units, unit_catalogue, unit_room
from mobilink.simulate import simulate_tables


def test_make_units_merges_within_d_and_keeps_edge_lengths():
    f = pd.DataFrame(
        {
            "molecule_id": ["m1", "m1", "m1", "m1", "m2"],
            "start": [100, 1_200, 1_900, 5_000, 0],
            "end": [1_000, 1_800, 3_000, 5_500, 700],
        }
    )
    u = make_units(f, d=500)
    assert len(u) == 3
    a = u.iloc[0]
    assert (a.start, a.end, a.n_calls, a.g1, a.g2) == (100, 3_000, 3, 900, 1_100)
    assert (u.iloc[1].start, u.iloc[1].n_calls) == (5_000, 1)
    assert u.iloc[2].molecule_id == "m2"
    assert len(make_units(f, d=50)) == 5


def test_unit_room_reduces_to_gene_level_room_for_single_genes():
    rng = np.random.default_rng(4)
    L = rng.integers(2_000, 40_000, 300)
    ids = pd.Index([f"r{i}" for i in range(L.size)])
    rows = []
    for i in range(0, L.size, 3):  # B calls far apart (> d) on every third molecule
        for st in range(200, int(L[i]) - 1_500, 6_000):
            rows.append((ids[i], st, st + int(rng.integers(500, 1_400))))
    fb = pd.DataFrame(rows, columns=["molecule_id", "start", "end"])
    a = rng.integers(400, 1_500, 50).astype(float)
    cat = Catalogue(s=a, g1=a, g2=a, w=np.full(a.size, 1 / a.size), complete=True)
    x_u = unit_room(L, ids, make_units(fb, d=1_000), cat, 0.8)
    x_g = excluded_volume_exposure(L, ids, fb, a, 0.8)
    assert np.allclose(x_u, x_g, rtol=1e-12, atol=1e-6)


def test_catalogue_uses_complete_units_with_inverse_probability_weights():
    lengths = pd.DataFrame({"molecule_id": ["a", "b"], "length": [20_000, 20_000]})
    u = pd.DataFrame(
        {
            "molecule_id": ["a", "a", "b"],
            "start": [100, 5_000, 2_000],
            "end": [900, 7_000, 12_000],
            "n_calls": [1, 2, 4],
            "g1": [800, 900, 1_000],
            "g2": [800, 1_000, 900],
        }
    )
    c = unit_catalogue(u, lengths, d=1_000)
    assert c.complete and c.s.tolist() == [2_000, 10_000]  # the unit at 100 bp from the start is not complete
    room = np.array([2 * (20_000 - 2_000 - 2_000), 2 * (20_000 - 10_000 - 2_000)], float)
    assert np.allclose(c.w, (1 / room) / (1 / room).sum())


def test_unit_nulls_are_near_one_for_clustered_b_either_way():
    rng = np.random.default_rng(8)
    reads = rng.lognormal(np.log(9_000), 0.7, 20_000).astype(int) + 500
    oe = {"A": [], "B": []}
    for seed in range(3):
        s = simulate_tables(
            reads,
            rng.integers(600, 1_300, 300),
            rng.integers(600, 1_500, 300),
            seed=100 + seed,
            n_genomes=100,
            total_bases=2e8,
            n_a=24.0,
            islands=5,
            singletons=12,
            clustered=True,
            linked=0.0,
            min_frac=0.8,
            n_truth_anchors=1000,
        )
        f, ln = s.features, s.lengths
        for fx in ("A", "B"):
            u = cooccurrence_test_units(
                ln,
                f[f["class"] == "ARG"],
                f[f["class"] == "MGE"],
                strata=ln.genome.to_numpy(),
                fixed=fx,
                tost_bound=None,
            )
            oe[fx].append(u.n2.observed / u.n2.expected)
    for fx in ("A", "B"):
        assert 0.85 < np.mean(oe[fx]) < 1.15, (fx, oe[fx])
