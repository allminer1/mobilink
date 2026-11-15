import numpy as np
import pandas as pd

from linkstat.features import covered_bases
from linkstat.nulls.analytic import cooccurrence_test
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
            total_bases=4e8,
            n_a=24.0,
            islands=5,
            singletons=12,
            linked=0.0,
            min_frac=0.8,
            n_truth_anchors=1_000,
        )
        f, ids, ln = sim.features, sim.lengths.molecule_id, sim.lengths.length.to_numpy()
        fa, fb = f[f["class"] == "ARG"], f[f["class"] == "MGE"]
        ha, hb = ids.isin(fa.molecule_id).to_numpy(), ids.isin(fb.molecule_id).to_numpy()
        la, lb = (fa.end - fa.start).to_numpy(float), (fb.end - fb.start).to_numpy(float)
        occ = covered_bases(fb, ids.to_numpy())
        oe["N0"].append(cooccurrence_test(ln, ha, hb, la, lb, tost_bound=None).oe)
        r = cooccurrence_test(ln, ha, hb, la, lb, tost_bound=None, occupied_b=occ)
        assert r.excluded_volume
        oe["N0ev"].append(r.oe)
        oe["N2"].append(binned_test(ln, ha, hb).oe)
        oe["N2ev"].append(binned_test(ln, ha, hb, bin_by=np.clip(ln - occ, 0, None)).oe)
    m = {k: np.mean(v) for k, v in oe.items()}
    assert m["N0"] < 0.88 and m["N2"] < 0.92  # length-only nulls: biased low by excluded volume
    assert abs(m["N0ev"] - 1) < 0.07 and abs(m["N2ev"] - 1) < 0.1
