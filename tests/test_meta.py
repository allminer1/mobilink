import numpy as np
import pandas as pd
from click.testing import CliRunner
from scipy import stats

from mobilink.cli import main
from mobilink.meta import contrasts, group_estimate, holm, jackknife_se
from mobilink.simulate import simulate_tables
from mobilink.stats.nettest import net_linked_test, net_w_estimate

G = np.random.default_rng(0)
READS = G.lognormal(np.log(6000), 0.7, 5000).astype(int) + 500
A_LEN = np.maximum(150, G.lognormal(np.log(900), 0.35, 2000)).astype(int)
B_LEN = np.maximum(150, G.lognormal(np.log(1000), 0.5, 2000)).astype(int)


def toy(linked, seed):
    return simulate_tables(READS, A_LEN, B_LEN, seed=seed, n_genomes=20, total_bases=2 * 20 * 4e6, islands=5,
                           singletons=12, clustered=True, linked=linked, min_frac=0.8, n_truth_anchors=5_000,
                           genome_length=(4_000_000, 4_000_000), density_sigma=1.0, truth_stratified=True,
                           positions=True, n_a=24.0)


def run_dict(name, linked, seed):
    s = toy(linked, seed)
    return dict(run=name, lengths=s.lengths, features=s.features, strata=s.lengths.genome.to_numpy())


def test_estimate_only_path_equals_the_test():
    s = toy(0.1, 21)
    ln, f, st = s.lengths, s.features, s.lengths.genome.to_numpy()
    fa, fb = f[f["class"] == "ARG"], f[f["class"] == "MGE"]
    for strata in (st, None):
        r = net_linked_test(ln, fa, fb, strata=strata, loci="positions", n_mc=99)
        e = net_w_estimate(ln, fa, fb, strata)
        assert e["net_w"] == r["net_w"] and e["units"] == r["units"] and e["sides_w"] > 0


def test_group_estimate_jackknife_and_contrasts():
    hi = [run_dict(f"h{i}", 0.25, 30 + i) for i in range(4)]
    lo = [run_dict(f"l{i}", 0.0, 40 + i) for i in range(3)]
    row_hi, rep_hi = group_estimate(hi, ["p1", "p1", "p2", "p3"], "ARG", "MGE")
    row_lo, rep_lo = group_estimate(lo, ["q1", "q2", "q3"], "ARG", "MGE")
    assert row_hi["clusters"] == 3 and list(rep_hi.left_out) == ["none", "p1", "p2", "p3"]
    assert rep_hi.runs.tolist() == [4, 2, 3, 3]
    reps = rep_hi.net_w.to_numpy()[1:]
    se = np.sqrt(2 / 3 * np.sum((reps - reps.mean()) ** 2))
    assert np.isclose(row_hi["se_jack"], se) and row_hi["df"] == 2
    assert np.isclose(row_hi["hi95"], row_hi["net_w"] + stats.t.ppf(0.975, 2) * se)
    assert 1 < row_hi["clusters_effective"] <= 3
    assert row_hi["net_w"] > 0.05 and abs(row_lo["net_w"]) < 0.03
    t = pd.DataFrame([dict(group="hi", **row_hi), dict(group="lo", **row_lo)])
    c = contrasts(t).iloc[0]
    assert c.a == "hi" and c.b == "lo" and np.isclose(c["diff"], row_hi["net_w"] - row_lo["net_w"])
    assert np.isclose(c.se, np.hypot(row_hi["se_jack"], row_lo["se_jack"])) and c.p_two_sided_holm == c.p_two_sided
    assert jackknife_se(0.1, [0.1])[0] != jackknife_se(0.1, [0.1])[0]  # NaN with one replicate
    assert np.allclose(holm([0.01, 0.04, 0.03]), [0.03, 0.06, 0.06])


def test_cli_meta(tmp_path):
    rows = []
    for i, (grp, linked) in enumerate([("gut", 0.25), ("gut", 0.25), ("gut", 0.25), ("soil", 0.0), ("soil", 0.0)]):
        s = toy(linked, 50 + i)
        s.lengths.to_csv(tmp_path / f"r{i}.lengths.tsv", sep="\t", index=False)
        s.features.to_csv(tmp_path / f"r{i}.features.tsv", sep="\t", index=False)
        rows.append(dict(run=f"r{i}", env=grp, project=f"{grp}{i}", lengths=f"r{i}.lengths.tsv",
                         features=f"r{i}.features.tsv"))
    pd.DataFrame(rows).to_csv(tmp_path / "manifest.tsv", sep="\t", index=False)
    out = tmp_path / "survey"
    man = str(tmp_path / "manifest.tsv")
    r = CliRunner().invoke(main, ["meta", "--manifest", man, "--pair", "ARG:MGE", "--group-col", "env", "--cluster-col",
                                  "project", "--strata-col", "genome", "--out", str(out)])
    assert r.exit_code == 0, r.output
    g = pd.read_csv(f"{out}.groups.tsv", sep="\t").set_index("env")
    j = pd.read_csv(f"{out}.jackknife.tsv", sep="\t")
    c = pd.read_csv(f"{out}.contrasts.tsv", sep="\t")
    assert g.loc["gut", "clusters"] == 3 and g.loc["soil", "runs"] == 2
    assert g.loc["gut", "net_w"] > g.loc["soil", "net_w"]
    assert len(j) == (1 + 3) + (1 + 2) and len(c) == 1 and c.a[0] == "gut"
    r = CliRunner().invoke(main, ["meta", "--manifest", man, "--pair", "ARG:MGE", "--out", str(out)])
    assert r.exit_code != 0 and "lacks the column" in r.output
