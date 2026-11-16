import json

import numpy as np
import pandas as pd
from click.testing import CliRunner

from linkstat.cli import main
from linkstat.pairs import net_test_pairs
from linkstat.simulate import simulate_tables
from linkstat.stats.nettest import net_linked_test

G = np.random.default_rng(0)
READS = G.lognormal(np.log(6000), 0.7, 8000).astype(int) + 500
A_LEN = np.maximum(150, G.lognormal(np.log(900), 0.35, 2000)).astype(int)
B_LEN = np.maximum(150, G.lognormal(np.log(1000), 0.5, 2000)).astype(int)


def toy(linked, seed=11):
    return simulate_tables(READS, A_LEN, B_LEN, seed=seed, n_genomes=30, total_bases=2 * 30 * 4e6, islands=5,
                           singletons=12, clustered=True, linked=linked, min_frac=0.8, n_truth_anchors=20_000,
                           genome_length=(4_000_000, 4_000_000), density_sigma=1.0, truth_stratified=True,
                           positions=True, n_a=24.0)


def test_net_test_pairs_matches_the_registered_function_and_adds_q_values():
    s = toy(0.1)
    ln, f = s.lengths, s.features
    st = ln.genome.to_numpy()
    t = net_test_pairs(ln, f, [("ARG", "MGE"), ("ARG", "XXX")], strata=st, loci="positions", n_mc=499)
    r = net_linked_test(ln, f[f["class"] == "ARG"], f[f["class"] == "MGE"], strata=st, loci="positions", n_mc=499)
    row = t.iloc[0]
    for k in ("net_w", "net_w_lo", "net_w_p", "net_mol", "net_mol_p"):
        assert row[k] == r[k] or (np.isnan(row[k]) and np.isnan(r[k]))
    assert row.net_w > 0 and row.net_w_p < 0.01 and row.net_w_q <= 0.01 + 1e-12
    assert t.iloc[1]["flags"] == "class empty" and np.isnan(t.iloc[1]["net_w_q"])
    assert "oe_n2u" in t.columns and row.oe_n2u > 1


def test_exclude_drops_molecules_before_the_test():
    s = toy(0.0, seed=12)
    ln, f = s.lengths, s.features
    drop = set(ln.molecule_id.iloc[:500])
    t = net_test_pairs(ln, f, [("ARG", "MGE")], strata=ln.genome.to_numpy(), loci="positions", exclude=drop,
                       n_mc=199, oe=False)
    assert t.iloc[0].excluded_molecules == 500 and t.iloc[0].molecules == len(ln) - 500


def test_cli_test_with_strata_column_and_exclude(tmp_path):
    s = toy(0.0, seed=13)
    s.lengths.to_csv(tmp_path / "l.tsv", sep="\t", index=False)
    s.features.to_csv(tmp_path / "f.tsv", sep="\t", index=False)
    pd.DataFrame({"read_id": s.lengths.molecule_id.iloc[:10]}).to_csv(tmp_path / "x.tsv", sep="\t", index=False)
    out = tmp_path / "res.tsv"
    args = ["test", "--lengths", str(tmp_path / "l.tsv"), "--features", str(tmp_path / "f.tsv"), "--pairs", "ARG:MGE",
            "--strata-col", "genome", "--exclude", str(tmp_path / "x.tsv"), "--n-mc", "199", "--out", str(out)]
    r = CliRunner().invoke(main, args)
    assert r.exit_code == 0, r.output
    t = pd.read_csv(out, sep="\t")
    assert len(t) == 1 and t.excluded_molecules.iloc[0] == 10 and 0 <= t.net_w_p.iloc[0] <= 1
    assert "no locus clumps" in t["flags"].iloc[0]
    meta = json.loads((tmp_path / "res.tsv.meta.json").read_text())
    assert meta["strata"] == "column genome" and meta["excluded_molecules"] == 10
    bad = CliRunner().invoke(main, args[:-2] + ["--kraken", str(tmp_path / "x.tsv"), "--out", str(out)])
    assert bad.exit_code != 0


def test_same_orf_default_is_the_registered_rule_and_none_turns_it_off():
    s = toy(0.0, seed=14)
    ln, f = s.lengths, s.features
    fa, fb = f[f["class"] == "ARG"], f[f["class"] == "MGE"]
    dual = fa.iloc[:5].assign(**{"class": "MGE"})  # five A calls also called as B (one ORF, two databases)
    fb2 = pd.concat([fb, dual], ignore_index=True)
    st = ln.genome.to_numpy()
    r0 = net_linked_test(ln, fa, fb2, strata=st, loci="positions", n_mc=99)
    r1 = net_linked_test(ln, fa, fb2, strata=st, loci="positions", n_mc=99, same_orf=0.5)
    r2 = net_linked_test(ln, fa, fb2, strata=st, loci="positions", n_mc=99, same_orf=None)
    assert r0["net_w"] == r1["net_w"] and r0["n_dual"] == 5 and r2["n_dual"] == 0
