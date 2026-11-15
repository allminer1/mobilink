import numpy as np
import pandas as pd

from linkstat.loci import accepted_overlaps, clusters_from_positions, locus_clusters
from linkstat.nulls.binned import binned_test
from linkstat.nulls.robust import cluster_robust_test


def test_mge_guard_rejects_overlaps_inside_an_element():
    paf = pd.DataFrame(
        {
            "q": ["r1", "r1", "r3"],
            "qlen": [20_000] * 3,
            "qs": [1_000, 10_000, 0],
            "qe": [3_500, 15_000, 4_000],
            "strand": ["+"] * 3,
            "t": ["r2", "r3", "r4"],
            "tlen": [20_000] * 3,
            "ts": [5_000, 0, 6_000],
            "te": [7_500, 5_000, 10_000],
            "nmatch": [2_000, 4_500, 3_600],
            "block": [2_500, 5_000, 4_000],
            "mapq": [0] * 3,
        }
    )
    mge = pd.DataFrame(
        {"molecule_id": ["r1", "r2", "r3"], "start": [800, 4_900, 0], "end": [3_600, 7_600, 1_000]}
    )  # r1-r2 overlap lies inside an IS element on both reads
    acc = accepted_overlaps(paf, mge, min_free=1_000)
    assert list(zip(acc.q, acc.t, strict=True)) == [("r1", "r3"), ("r3", "r4")]
    lab = locus_clusters(["r1", "r2", "r3", "r4", "r5"], paf, mge)
    assert lab[0] == lab[2] == lab[3] and len({lab[0], lab[1], lab[4]}) == 3


def test_clusters_from_positions():
    rep = [0, 0, 0, 0, 1]
    start = [0, 7_000, 28_000, 20_500, 7_000]
    length = [10_000, 10_000, 5_000, 10_000, 10_000]
    lab = clusters_from_positions(rep, start, length, min_block=2_000)
    assert lab[0] == lab[1] and lab[2] == lab[3] and len(set(lab)) == 3


def test_clusters_from_positions_are_not_broken_by_short_reads():
    # a long chain with short reads in between: the short reads are singletons, the chain stays one cluster
    rep = [0] * 7
    start = [0, 3_000, 3_500, 6_000, 6_200, 9_000, 30_000]
    length = [8_000, 500, 8_000, 300, 8_000, 8_000, 5_000]
    lab = clusters_from_positions(rep, start, length, min_block=2_000)
    assert lab[0] == lab[2] == lab[4] == lab[5]
    assert len({lab[0], lab[1], lab[3], lab[6]}) == 4


def _null_data(rng, n=40_000):
    L = rng.integers(2_000, 30_000, n)
    pa = 1 - np.exp(-L / 150_000)
    pb = 1 - np.exp(-L / 40_000)
    return L, rng.random(n) < pa, rng.random(n) < pb


def test_design_effect_is_about_one_without_and_k_with_k_copies():
    rng = np.random.default_rng(5)
    L, a, b = _null_data(rng)
    r1 = cluster_robust_test(a, b, L, np.arange(L.size))
    assert 0.8 < r1.deff < 1.25
    k = 5  # every molecule read k times: one cluster per locus
    Lk, ak, bk = np.repeat(L, k), np.repeat(a, k), np.repeat(b, k)
    rk = cluster_robust_test(ak, bk, Lk, np.repeat(np.arange(L.size), k))
    assert 0.8 * k < rk.deff < 1.25 * k
    ex = binned_test(Lk, ak, bk, bin_by=Lk)
    rk2 = cluster_robust_test(ak, bk, Lk, np.repeat(np.arange(L.size), k), exact=ex)
    assert rk2.p_two_sided >= ex.p_two_sided and rk2.p_two_sided >= rk2.p_two_cr
