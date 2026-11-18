import numpy as np
import pytest

from mobilink.power import (
    calls_density,
    depth_table,
    detectable_oe,
    netlinked_halfwidth,
    power_at,
    required_bases,
    sides_reaching,
)

RNG = np.random.default_rng(0)
L = np.exp(RNG.normal(np.log(8_000), 0.7, 5_000)).astype(int) + 300
GENE = np.array([1_000.0])


def test_power_at_detectable_oe_is_the_target_power():
    d = detectable_oe(L, GENE, GENE, 1 / 30_000, 1 / 15_000, 50_000, power=0.8, min_frac=1.0)
    assert d["oe_enrichment"] > 1 > d["oe_depletion"]
    for r in (d["oe_enrichment"], d["oe_depletion"]):
        assert power_at(L, GENE, GENE, 1 / 30_000, 1 / 15_000, 50_000, r, min_frac=1.0) == pytest.approx(0.8, abs=0.01)
    # more molecules, smaller detectable effect
    d2 = detectable_oe(L, GENE, GENE, 1 / 30_000, 1 / 15_000, 500_000, power=0.8, min_frac=1.0)
    assert d2["oe_enrichment"] < d["oe_enrichment"]


def test_sides_reaching_matches_monte_carlo_placement():
    rng = np.random.default_rng(1)
    lam, g, w = 1 / 20_000, 1_000, 2_000
    count = 0
    for Li in L:
        n = rng.poisson(lam * max(Li - g, 0))
        s = rng.random(n) * (Li - g)
        count += int(np.sum(Li - (s + g) >= w) + np.sum(s >= w))
    assert sides_reaching(L, np.array([g]), lam, w) == pytest.approx(count, rel=0.05)


def test_halfwidth_and_required_depth_are_consistent():
    assert netlinked_halfwidth(100, 0.1, 1 / 5_000, 2_000) > netlinked_halfwidth(400, 0.1, 1 / 5_000, 2_000)
    lam_a = calls_density(3_000, L, GENE)
    g = required_bases(L, GENE, lam_a, 1 / 5_000, 2_000, 0.05)
    t = depth_table(L, GENE, GENE, 1 / 30_000, 1 / 15_000, lam_a, 1 / 5_000, [g], min_frac=1.0)
    assert t.halfwidth.iloc[0] == pytest.approx(0.05, rel=1e-4)


def test_at_risk_and_power_cli(tmp_path):
    import pandas as pd
    from click.testing import CliRunner

    from mobilink.cli import main
    from mobilink.power import at_risk, chance_and_linked
    from mobilink.simulate import simulate_tables

    f_r, f_a = chance_and_linked(0.1, 1 / 5_000, 2_000)
    assert at_risk(1_000, 0.1, 1 / 5_000, 2_000) == pytest.approx(1_000 * (1 - f_a))
    sim = simulate_tables(
        L,
        RNG.integers(600, 1_500, 200),
        RNG.integers(400, 2_000, 200),
        seed=2,
        n_genomes=20,
        total_bases=1e8,
        n_truth_anchors=20_000,
    )
    sim.lengths.to_csv(tmp_path / "l.tsv", sep="\t", index=False)
    sim.features.to_csv(tmp_path / "f.tsv", sep="\t", index=False)
    r = CliRunner().invoke(
        main,
        [
            "power",
            "--lengths",
            str(tmp_path / "l.tsv"),
            "--features",
            str(tmp_path / "f.tsv"),
            "--a",
            "ARG",
            "--b",
            "MGE",
            "--bases",
            "1e8",
            "--bases",
            "1e9",
        ],
    )
    assert r.exit_code == 0, r.output
    assert "half-width" in r.output and pd is not None
