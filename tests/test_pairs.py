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


def test_no_events_gives_no_estimate_instead_of_zero():
    s = toy(0.0, seed=15)
    ln, f = s.lengths.copy(), s.features
    hosts = ln.genome.to_numpy()
    fa = f[f["class"] == "ARG"]
    a_hosts = set(hosts[ln.molecule_id.isin(fa.molecule_id).to_numpy()])
    # one B call on a molecule of a host that carries no A: no B is ever near an A side or an anchor
    lonely = ln[~np.isin(hosts, list(a_hosts))].molecule_id
    if lonely.empty:  # every host carries A: put the B on a molecule in a new stratum
        ln.loc[ln.index[-1], "genome"] = "solo"
        lonely = ln.molecule_id.iloc[[-1]]
        fa = fa[fa.molecule_id != lonely.iloc[0]]
    fb = pd.DataFrame({"molecule_id": [lonely.iloc[0]], "start": [10], "end": [800], "class": ["MRG"]})
    t = net_test_pairs(ln, pd.concat([fa, fb]), [("ARG", "MRG")], strata=ln.genome.to_numpy(), loci="positions",
                       n_mc=99, oe=False).iloc[0]
    assert np.isnan(t.net_w) and np.isnan(t.net_w_lo) and "net_w: no events" in t["flags"]


def test_mol_ratio_is_reported_with_a_90_interval():
    s = toy(0.1, seed=16)
    ln, f = s.lengths, s.features
    r = net_linked_test(ln, f[f["class"] == "ARG"], f[f["class"] == "MGE"], strata=ln.genome.to_numpy(),
                        loci="positions", n_mc=99)
    assert np.isclose(r["mol_ratio"], r["mol_share_A"] / r["mol_share_R"])
    assert np.isclose(r["net_mol"], 1 - (1 - r["mol_share_A"]) / (1 - r["mol_share_R"]))
    assert r["mol_ratio_lo90"] < r["mol_ratio"] < r["mol_ratio_hi90"] and r["mol_ratio"] > 1


def test_pool_runs_keeps_runs_apart_and_flags_unresolved():
    from linkstat.pairs import pool_runs

    runs = []
    for k, seed in enumerate((17, 18)):
        s = toy(0.1, seed=seed)
        st = s.lengths.genome.astype(str).to_numpy().astype(object)
        st[:200] = "unresolved"
        runs.append(dict(run=f"run{k}", lengths=s.lengths, features=s.features, strata=st))
    ln, f, st, paf = pool_runs(runs)
    assert paf is None and ln.molecule_id.is_unique and len(ln) == sum(len(r["lengths"]) for r in runs)
    assert set(x.split("|")[0] for x in st) == {"run0", "run1"}
    t = net_test_pairs(ln, f, [("ARG", "MGE")], strata=st, n_mc=199, oe=False).iloc[0]
    assert t.net_w > 0 and t.net_w_p < 0.05 and "a_molecules_unresolved" in t.index
    one = net_test_pairs(runs[0]["lengths"], runs[0]["features"], [("ARG", "MGE")], strata=runs[0]["strata"],
                         n_mc=199, oe=False).iloc[0]
    assert t.units > one.units


def test_read_paf_keep_filter(tmp_path):
    from linkstat.loci import PAF_COLS, read_paf

    rows = [("a", 5000, 0, 3000, "+", "b", 6000, 100, 3100, 2900, 3000, 60),
            ("a", 5000, 0, 3000, "+", "c", 6000, 100, 3100, 2900, 3000, 60),
            ("a", 5000, 0, 500, "+", "b", 6000, 100, 600, 480, 500, 60)]
    pd.DataFrame(rows, columns=PAF_COLS).to_csv(tmp_path / "x.paf", sep="\t", header=False, index=False)
    assert len(read_paf(tmp_path / "x.paf", min_block=1000)) == 2
    assert len(read_paf(tmp_path / "x.paf", min_block=1000, keep={"a", "b"})) == 1
    (tmp_path / "empty.paf").write_text("")
    assert len(read_paf(tmp_path / "empty.paf")) == 0


def test_molecule_clumps_keep_the_estimates_and_widen_the_test(tmp_path):
    s = toy(0.1, seed=13)
    ln, f = s.lengths, s.features
    st = ln.genome.to_numpy()
    none = net_test_pairs(ln, f, [("ARG", "MGE")], strata=st, n_mc=499, oe=False).iloc[0]
    mol = net_test_pairs(ln, f, [("ARG", "MGE")], strata=st, loci="molecule", n_mc=499, oe=False).iloc[0]
    for k in ("net_w", "net_mol", "net_w_O", "net_w_E", "net_mol_O", "net_mol_E"):
        assert mol[k] == none[k]
    assert mol.unit_loci == f[f["class"] == "ARG"].molecule_id.nunique() <= none.unit_loci
    assert mol.net_w_V >= none.net_w_V and mol.net_mol_V >= none.net_mol_V
    assert "molecule clumps" in mol["flags"] and "locus conflicts" not in mol["flags"]
    assert "no locus clumps" in none["flags"] and "no locus clumps" not in mol["flags"]
    ln.to_csv(tmp_path / "l.tsv", sep="\t", index=False)
    f.to_csv(tmp_path / "f.tsv", sep="\t", index=False)
    args = ["test", "--lengths", str(tmp_path / "l.tsv"), "--features", str(tmp_path / "f.tsv"), "--pairs", "ARG:MGE",
            "--strata-col", "genome", "--n-mc", "499", "--no-oe", "--out", str(tmp_path / "res.tsv")]
    r = CliRunner().invoke(main, args + ["--molecule-clumps"])
    assert r.exit_code == 0, r.output
    row = pd.read_csv(tmp_path / "res.tsv", sep="\t").iloc[0]
    assert np.isclose(row.net_w, mol.net_w) and np.isclose(row.net_w_p, mol.net_w_p)
    assert json.load(open(tmp_path / "res.tsv.meta.json"))["molecule_clumps"] is True
    empty = tmp_path / "ava.paf"
    empty.write_text("")
    bad = CliRunner().invoke(main, args + ["--molecule-clumps", "--loci", str(empty)])
    assert bad.exit_code != 0 and "either --loci" in bad.output
