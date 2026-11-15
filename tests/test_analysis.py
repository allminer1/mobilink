import json

import numpy as np
import pandas as pd
from click.testing import CliRunner

from linkstat.analysis import analyse_pairs
from linkstat.cli import main
from linkstat.simulate import simulate_tables

RNG = np.random.default_rng(3)
READS = RNG.lognormal(np.log(9_000), 0.6, 10_000).astype(int) + 500
A_LEN = RNG.integers(600, 1_500, 300)
B_LEN = RNG.integers(400, 2_000, 300)


def toy(linked=0.3, seed=1):
    return simulate_tables(
        READS,
        A_LEN,
        B_LEN,
        seed=seed,
        n_genomes=40,
        total_bases=3e8,
        linked=linked,
        min_frac=0.8,
        n_truth_anchors=50_000,
    )


def test_analyse_pairs_recovers_planted_linkage_and_reports_all_nulls():
    sim = toy()
    lengths = sim.lengths.assign(sample=np.where(np.arange(len(sim.lengths)) % 2 == 0, "s1", "s2"))
    feats = pd.concat(
        [
            sim.features,
            pd.DataFrame(
                {"molecule_id": ["sim0"], "start": [0], "end": [100], "class": ["MRG"], "scov": [100.0], "locus": [-1]}
            ),
        ]
    )
    t = analyse_pairs(lengths, feats, [("ARG", "MGE"), ("ARG", "XXX")], strata_col="sample", n_boot=50)
    r = t.iloc[0]
    truth = float(sim.truth.query("mode == 'end' and w == 2000").net.iloc[0])
    assert r.intact_lo - 0.03 < truth < r.intact_hi + 0.03
    assert r.oe_ratio > 1.5 and r.p < 1e-6 and r.null_model == "N2-ev"
    for c in ("registered_oe", "n0ev_oe", "mh_or", "p_mid", "p_tost", "naive_share", "q", "flags"):
        assert c in t.columns
    assert "pseudo-replication" in r["flags"]
    assert t.iloc[1]["flags"] == "class empty"


def test_cli_test_writes_results_and_meta(tmp_path):
    sim = toy(linked=0.0, seed=2)
    sim.lengths.to_csv(tmp_path / "l.tsv", sep="\t", index=False)
    sim.features.to_csv(tmp_path / "f.tsv", sep="\t", index=False)
    out = tmp_path / "res.tsv"
    r = CliRunner().invoke(
        main,
        [
            "test",
            "--lengths",
            str(tmp_path / "l.tsv"),
            "--features",
            str(tmp_path / "f.tsv"),
            "--pairs",
            "ARG:MGE",
            "--boot",
            "20",
            "--out",
            str(out),
        ],
    )
    assert r.exit_code == 0, r.output
    t = pd.read_csv(out, sep="\t")
    assert len(t) == 1 and abs(t.intact_net.iloc[0]) < 0.06
    meta = json.loads((tmp_path / "res.tsv.meta.json").read_text())
    assert meta["pairs"] == ["ARG:MGE"] and len(meta["lengths_sha256"]) == 64
    bad = CliRunner().invoke(
        main,
        [
            "test",
            "--lengths",
            str(tmp_path / "l.tsv"),
            "--features",
            str(tmp_path / "f.tsv"),
            "--pairs",
            "ARGMGE",
            "--out",
            str(out),
        ],
    )
    assert bad.exit_code != 0
