import numpy as np
import pandas as pd
from click.testing import CliRunner

from linkstat.cli import main
from linkstat.io.tables import validate_features
from linkstat.simulate import simulate_tables
from linkstat.stats.netlinked import length_robustness

RNG = np.random.default_rng(0)
READS = RNG.lognormal(np.log(9_000), 0.6, 5_000).astype(int) + 500
A_LEN = RNG.integers(600, 1_500, 300)
B_LEN = RNG.integers(400, 2_000, 300)


def truth_net(sim, w=2_000):
    return float(sim.truth.query("mode == 'end' and w == @w").net.iloc[0])


def test_tables_are_valid_and_carry_loci():
    sim = simulate_tables(
        READS, A_LEN, B_LEN, seed=1, n_genomes=10, total_bases=2e7, min_frac=0.8, n_truth_anchors=20_000
    )
    validate_features(sim.features, sim.lengths)
    assert set(sim.features["class"]) == {"ARG", "MGE"}
    assert sim.features.locus.notna().all() and sim.features.locus.nunique() < len(sim.features)


def test_truth_is_zero_without_linkage_and_rises_with_it():
    nets = [
        truth_net(
            simulate_tables(
                READS,
                A_LEN,
                B_LEN,
                seed=2,
                n_genomes=20,
                total_bases=1e6,
                linked=p,
                clustered=False,
                n_truth_anchors=100_000,
            )
        )
        for p in (0.0, 0.3, 0.6)
    ]
    assert abs(nets[0]) < 0.02 and nets[0] < nets[1] < nets[2]


def test_estimator_recovers_truth_on_simulated_reads():
    sim = simulate_tables(
        READS, A_LEN, B_LEN, seed=3, n_genomes=40, total_bases=4e8, linked=0.4, min_frac=0.8, n_truth_anchors=100_000
    )
    f = sim.features
    row = length_robustness(
        sim.lengths, f[f["class"] == "ARG"], f[f["class"] == "MGE"], targets=(), n_boot=50, n_anchors=10_000
    )
    assert abs(row["intact_net"] - truth_net(sim)) < 0.04


def test_cli_writes_three_tables(tmp_path):
    p = tmp_path / "rl.tsv"
    pd.DataFrame({"read_id": [f"r{i}" for i in range(READS.size)], "length": READS}).to_csv(p, sep="\t", index=False)
    r = CliRunner().invoke(
        main,
        [
            "simulate",
            "--read-lengths",
            str(p),
            "--n-genomes",
            "8",
            "--total-bases",
            "5e6",
            "--out-prefix",
            str(tmp_path / "sim"),
        ],
    )
    assert r.exit_code == 0, r.output
    assert (tmp_path / "sim.features.tsv.gz").exists() and (tmp_path / "sim.truth.tsv").exists()
    assert "true net linked" in r.output
