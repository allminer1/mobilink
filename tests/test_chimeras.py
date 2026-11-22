import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from mobilink.chimeras import LSK_CORE, adapter_hits, chimera_flags, overlap_breakpoints, overlap_breakpoints_stream
from mobilink.cli import main
from mobilink.loci import PAF_COLS


def paf_rows():
    # chimeric read R (10 kb, junction at 5 kb): its left half matches the middle of P1 and P2, its right half the
    # middle of P3; all partners continue past the junction, so the overlaps end inside both reads there.
    # clean read C: dovetail overlaps that end where the partner ends, and one read spanning 3-7 kb.
    r = [
        ("R", 10000, 0, 5000, "+", "P1", 12000, 3000, 8000, 4800, 5000, 60),
        ("R", 10000, 100, 5000, "+", "P2", 12000, 4000, 8900, 4700, 4900, 60),
        ("R", 10000, 5000, 10000, "+", "P3", 12000, 3000, 8000, 4800, 5000, 60),
        ("C", 10000, 0, 6000, "+", "Q1", 9000, 3000, 9000, 5800, 6000, 60),
        ("C", 10000, 4000, 10000, "+", "Q2", 9000, 0, 6000, 5800, 6000, 60),
        ("C", 10000, 2000, 9000, "-", "Q3", 7000, 0, 7000, 6800, 7000, 60),
    ]
    return pd.DataFrame(r, columns=PAF_COLS)


def test_overlap_breakpoint_on_a_chimera_and_none_on_a_clean_read():
    assessed, bp = overlap_breakpoints(paf_rows())
    assert {"R", "C", "P1", "P2", "P3", "Q1", "Q2", "Q3"} <= assessed
    r = bp[bp.read_id == "R"]
    assert len(r) == 1 and abs(r.x.iloc[0] - 5000) <= 200 and r.support.iloc[0] == 3
    assert (bp.read_id == "C").sum() == 0
    t = chimera_flags(paf=paf_rows())
    assert t.set_index("read_id").loc["R", "flag"] and not t.set_index("read_id").loc["C", "flag"]
    t1 = chimera_flags(paf=paf_rows().iloc[[0, 3, 4, 5]], min_support=2)  # one partner only: not enough support
    assert not t1.set_index("read_id").loc["R", "flag"]


def test_adapter_inside_a_read_is_found_but_not_at_its_ends():
    pytest.importorskip("edlib")
    rng = np.random.default_rng(1)
    def rnd(n):
        return "".join(rng.choice(list("ACGT"), n))
    reads = [("mid", rnd(3000) + LSK_CORE + rnd(3000)), ("end", LSK_CORE + rnd(6000)), ("none", rnd(6000))]
    h = adapter_hits(reads)
    assert set(h.read_id) == {"mid"} and (h.pattern == "lsk").all()


def test_cli_chimeras_table_feeds_test_exclude(tmp_path):
    p = paf_rows()
    p.to_csv(tmp_path / "ava.paf", sep="\t", header=False, index=False)
    out = tmp_path / "flags.tsv"
    r = CliRunner().invoke(main, ["chimeras", "--paf", str(tmp_path / "ava.paf"), "--out", str(out)])
    assert r.exit_code == 0, r.output
    f = pd.read_csv(out, sep="\t")
    assert f.loc[f.flag, "read_id"].tolist() == ["R"]
    from mobilink.simulate import simulate_tables

    g = np.random.default_rng(2)
    sim = simulate_tables(g.lognormal(np.log(6000), 0.6, 3000).astype(int) + 500, g.integers(600, 1500, 200),
                          g.integers(400, 2000, 200), seed=3, n_genomes=20, total_bases=5e7, linked=0.0,
                          min_frac=0.8, n_truth_anchors=5_000)
    ids = sim.lengths.molecule_id
    flags = pd.DataFrame({"read_id": ids.iloc[:30], "flag": [True] * 10 + [False] * 20})
    flags.to_csv(tmp_path / "flags2.tsv", sep="\t", index=False)
    sim.lengths.to_csv(tmp_path / "l.tsv", sep="\t", index=False)
    sim.features.to_csv(tmp_path / "f.tsv", sep="\t", index=False)
    r = CliRunner().invoke(main, ["test", "--lengths", str(tmp_path / "l.tsv"), "--features", str(tmp_path / "f.tsv"),
                                  "--pairs", "ARG:MGE", "--exclude", str(tmp_path / "flags2.tsv"), "--n-mc", "99",
                                  "--no-oe", "--out", str(tmp_path / "res.tsv")])
    assert r.exit_code == 0, r.output
    assert pd.read_csv(tmp_path / "res.tsv", sep="\t").excluded_molecules.iloc[0] == 10


def test_tiny_input_is_reported_not_raised():
    from mobilink.pairs import net_test_pairs

    reads = pd.DataFrame({"molecule_id": ["R", "C", "X"], "length": [10000, 10000, 8000]})
    feats = pd.DataFrame({"molecule_id": ["R", "C", "X", "X"], "start": [100, 100, 100, 3000],
                          "end": [900, 900, 900, 3900], "class": ["ARG", "ARG", "ARG", "MGE"]})
    t = net_test_pairs(reads, feats, [("ARG", "MGE")], n_mc=99, oe=False)
    assert len(t) == 1 and isinstance(t["flags"].iloc[0], str)


def random_paf(n_reads=400, n_overlaps=6000, seed=5):
    """Random overlaps among reads of fixed lengths; coordinates are valid but arbitrary (the two functions must agree
    on any table, sensible or not)."""
    g = np.random.default_rng(seed)
    length = g.integers(3000, 30000, n_reads)
    q, t = g.integers(0, n_reads, n_overlaps), g.integers(0, n_reads, n_overlaps)
    rows = []
    for a, b in zip(q, t, strict=True):
        block = int(g.integers(500, min(length[a], length[b])))
        qs, ts = int(g.integers(0, length[a] - block + 1)), int(g.integers(0, length[b] - block + 1))
        rows.append((f"r{a:04d}", int(length[a]), qs, qs + block, "+-"[int(g.integers(2))], f"r{b:04d}", int(length[b]), ts,
                     ts + block, int(block * g.uniform(0.1, 0.9)), block, 60))
    return pd.DataFrame(rows, columns=PAF_COLS)


@pytest.mark.parametrize("chunksize,buckets", [(10**6, 64), (700, 5), (97, 1)])
def test_streamed_paf_gives_the_same_breakpoints_as_the_table_in_memory(tmp_path, chunksize, buckets):
    for k, paf in enumerate((paf_rows(), random_paf())):
        path = tmp_path / f"ava{k}.paf.gz"
        paf.to_csv(path, sep="\t", header=False, index=False)
        a0, b0 = overlap_breakpoints(paf)
        a1, b1 = overlap_breakpoints_stream(path, chunksize=chunksize, buckets=buckets, tmp_dir=tmp_path)
        assert a0 == a1 and len(b0) > 0
        pd.testing.assert_frame_equal(b0.sort_values(["read_id", "x"]).reset_index(drop=True), b1, check_dtype=False)
        pd.testing.assert_frame_equal(chimera_flags(paf=paf), chimera_flags(paf=str(path), tmp_dir=tmp_path))
    assert not list(tmp_path.glob("mobilink_chimeras_*"))  # temporary files are removed


def test_streamed_paf_without_overlaps(tmp_path):
    (tmp_path / "empty.paf").write_text("")
    assessed, bp = overlap_breakpoints_stream(tmp_path / "empty.paf")
    assert assessed == set() and bp.empty and list(bp.columns) == ["read_id", "x", "support", "edge"]


@pytest.mark.parametrize("first", [True, False])
def test_spanning_overlap_that_starts_where_another_starts_still_vetoes(first):
    # P4 covers R from 0 to 9 kb and so spans the junction at 5 kb; P1 starts at the same position and ends at the
    # junction. The veto must not depend on which of the two rows comes first.
    span = pd.DataFrame([("R", 10000, 0, 9000, "+", "P4", 9000, 0, 9000, 8800, 9000, 60)], columns=PAF_COLS)
    paf = pd.concat([span, paf_rows()] if first else [paf_rows(), span], ignore_index=True)
    _, bp = overlap_breakpoints(paf)
    assert (bp.read_id == "R").sum() == 0
