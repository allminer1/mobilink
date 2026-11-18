import numpy as np
import pandas as pd
import pytest

from mobilink.fragment import fragment_partial, fragment_whole, n50, random_bounds, rate_for_n50
from mobilink.stats.netlinked import anchors, cut_mask, km, length_robustness, p_whole, prepare, sides


def test_km_matches_hand_computation():
    # times 1(event) 2(censored) 3(event) 4(event) 5(censored): S(3.5) = (1 - 1/5) * (1 - 1/3) = 0.5333
    t = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    e = np.array([1, 0, 1, 1, 0])
    assert km(t, e, 3.5, min_at_risk=1) == pytest.approx(1 - 0.8 * (2 / 3))
    # equal weights give the unweighted estimate; events at w count
    assert km(t, e, 3.5, np.full(5, 2.0), min_at_risk=1) == pytest.approx(km(t, e, 3.5, min_at_risk=1))
    assert km(t, e, 4.0, min_at_risk=1) == pytest.approx(1 - 0.8 * (2 / 3) * (1 / 2))
    # undefined if too few sides reach w
    assert np.isnan(km(t, e, 4.5, min_at_risk=2))


def test_km_events_before_censoring_at_ties():
    t = np.array([2.0, 2.0, 3.0])
    e = np.array([0, 1, 1])
    assert km(t, e, 2.0, min_at_risk=1) == pytest.approx(1 / 3)


def test_p_whole_brute_force():
    rng = np.random.default_rng(1)
    L = rng.integers(100, 20_000, 500)
    for g in (0, 50, 999, 5_000, 19_999, 25_000):
        assert p_whole(L, g) == pytest.approx(np.clip(L - g, 0, None).sum() / L.sum())


def test_sides_event_censoring_and_cut_calls():
    lengths = np.array([10_000])
    rid = np.array([0, 0])
    start = np.array([1_000, 2_500])
    end = np.array([2_000, 3_000])
    is_a = np.array([True, False])
    is_b = ~is_a
    i, t, e = sides(lengths, rid, start, end, is_a, is_b)
    # right side: far end of B (3,000) minus A end (2,000); left side: censored at A start
    assert sorted(zip(t, e, strict=True)) == [(1_000.0, 0), (1_000.0, 1)]
    # the same B marked cut censors the right side at the molecule end instead
    i, t, e = sides(lengths, rid, start, end, is_a, is_b, cut=np.array([False, True]))
    assert sorted(zip(t, e, strict=True)) == [(1_000.0, 0), (8_000.0, 0)]


def test_cut_mask_rules():
    lengths = np.array([5_000])
    rid = np.zeros(4, dtype=int)
    start = np.array([0, 100, 4_960, 2_000])
    end = np.array([900, 1_000, 4_990, 3_000])
    scov = np.array([99.0, 99.0, 90.0, 80.0])
    # touches the start; 100 bp from it (too far); 10 bp from the end with cover < 95; interior
    assert list(cut_mask(lengths, rid, start, end, scov)) == [True, False, True, False]


def test_fragmentation_helpers():
    assert n50([10, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]) == 10
    lengths = np.full(2_000, 20_000)
    rate = rate_for_n50(lengths, 3_000, 7)
    off, b = random_bounds(lengths, rate, 7)
    assert abs(n50(np.diff(b)) - 3_000) < 300
    assert set(off) <= set(b)  # every molecule start is a boundary
    # a feature split 90/10 by a breakpoint is kept (truncated) by fragment_partial and dropped by fragment_whole
    bounds = (np.array([0]), np.array([0, 1_090, 10_000]))
    one = (np.array([10_000]), np.array([0]), np.array([190]), np.array([1_190]), np.array(["A"]))
    pl, piece, s, e, c = fragment_partial(*one, bounds)
    assert list(piece) == [0] and list(s) == [190] and list(e) == [1_090]
    assert len(fragment_whole(*one, bounds)[1]) == 0


def test_anchors_fit_on_molecules():
    rng = np.random.default_rng(3)
    lengths = rng.integers(500, 30_000, 300)
    rid, st, en = anchors(lengths, rng.integers(300, 3_000, 100), 5_000, rng)
    assert len(rid) == 5_000 and (st >= 0).all() and (en <= lengths[rid]).all()


def simulated(n_reads=3_000, read_len=20_000, gene=900, lam_b=1 / 4_000, p_link=0.0, seed=5):
    """A genes uniform on reads; B genes as a Poisson process; with probability p_link a B is placed 100 bp after A."""
    rng = np.random.default_rng(seed)
    rows = []
    for r in range(n_reads):
        nb = rng.poisson(lam_b * read_len)
        for s in rng.integers(0, read_len - gene, nb):
            rows.append((f"r{r}", s, s + gene, "B"))
        s = int(rng.integers(0, read_len - 3 * gene))
        rows.append((f"r{r}", s, s + gene, "A"))
        if rng.random() < p_link:
            rows.append((f"r{r}", s + gene + 100, s + 2 * gene + 100, "B"))
    f = pd.DataFrame(rows, columns=["molecule_id", "start", "end", "class"]).assign(scov=100.0)
    ln = pd.DataFrame({"molecule_id": [f"r{r}" for r in range(n_reads)], "length": read_len})
    return ln, f[f["class"] == "A"].drop(columns="class"), f[f["class"] == "B"].drop(columns="class")


def test_net_linked_near_zero_without_linkage_and_positive_with_it():
    ln, fa, fb = simulated(p_link=0.0)
    row = length_robustness(ln, fa, fb, w=2_000, targets=(), n_boot=50, n_anchors=10_000)
    assert abs(row["intact_net"]) < 0.04
    assert row["intact_lo"] < 0 < row["intact_hi"] or abs(row["intact_net"]) < 0.02
    # the net linked fraction is per side: a B planted on the right side of 40% of A genes gives 0.4 / 2 = 0.2
    ln, fa, fb = simulated(p_link=0.4)
    row = length_robustness(ln, fa, fb, w=2_000, targets=(), n_boot=50, n_anchors=10_000)
    assert abs(row["intact_net"] - 0.2) < 0.04
    assert row["intact_lo"] < 0.2 < row["intact_hi"]


def test_length_robustness_keys_and_fragmented_versions():
    ln, fa, fb = simulated(n_reads=1_500, p_link=0.3)
    row = length_robustness(ln, fa, fb, w=2_000, targets=(5_000,), draws=2, n_boot=20, n_anchors=5_000)
    for k in ("intact_net", "intact_lo", "intact_hi", "r5000_net", "r5000_half", "r5000_sides_w", "n_dual"):
        assert k in row
    assert row["r5000_n50"] < row["intact_n50"]


def test_prepare_rejects_unknown_molecules():
    ln = pd.DataFrame({"molecule_id": ["r1"], "length": [5_000]})
    fa = pd.DataFrame({"molecule_id": ["r1"], "start": [0], "end": [900], "scov": [100.0]})
    fb = pd.DataFrame({"molecule_id": ["r2"], "start": [0], "end": [900], "scov": [100.0]})
    with pytest.raises(ValueError):
        prepare(ln, fa, fb, n_anchors=100)


def test_zero_event_upper_is_clopper_pearson():
    from scipy.stats import beta

    from mobilink.stats.netlinked import events_within, zero_event_upper

    for n in (1, 5, 30, 400):
        assert zero_event_upper(n) == pytest.approx(beta.ppf(0.975, 1, n))
    assert zero_event_upper(0) == 1.0
    assert events_within(np.array([1.0, 2.0, 3.0]), np.array([1, 0, 1]), 2.5) == 1
