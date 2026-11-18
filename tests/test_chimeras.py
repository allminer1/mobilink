import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from mobilink.chimeras import LSK_CORE, adapter_hits, chimera_flags, overlap_breakpoints
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
