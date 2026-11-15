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


def _p2_reads(seed=0):
    rng = np.random.default_rng(seed)
    return (
        rng.lognormal(np.log(10_000), 0.6, 20_000).astype(int) + 500,
        rng.integers(600, 1_300, 300),
        rng.integers(600, 1_500, 300),
    )


def test_a_islands_make_a_clustered():
    from linkstat.nulls.analytic import clustering_index

    reads, al, bl = _p2_reads()
    kw = dict(seed=9, n_genomes=60, total_bases=4e8, linked=0.0, min_frac=0.8, n_truth_anchors=1000)
    idx = {}
    for name, extra in (("flat", {}), ("islands", dict(n_a=6.0, a_islands=5.0))):
        s = simulate_tables(reads, al, bl, **kw, **extra)
        f, ln = s.features, s.lengths
        idx[name] = clustering_index(ln.length.to_numpy(), ln.molecule_id, f[f["class"] == "ARG"])
    assert idx["flat"] < 1.3 < 3 < idx["islands"], idx


def test_density_sigma_spreads_gene_counts_between_genomes():
    reads, al, bl = _p2_reads(1)
    cv = {}
    for sig in (0.0, 1.0):
        s = simulate_tables(
            reads,
            al,
            bl,
            seed=11,
            n_genomes=80,
            total_bases=8e8,
            linked=0.0,
            density_sigma=sig,
            genome_length=(4_000_000, 4_000_000),
            abundance_sigma=0.0,
            n_truth_anchors=1000,
        )
        f = s.features.merge(s.lengths[["molecule_id", "genome"]], on="molecule_id")
        per_g = f[f["class"] == "ARG"].groupby("genome").locus.nunique().reindex(range(80), fill_value=0)
        cv[sig] = per_g.std() / per_g.mean()
    assert cv[1.0] > 2 * cv[0.0], cv


def test_plasmids_are_replicons_labelled_with_their_host():
    reads, al, bl = _p2_reads(2)
    s = simulate_tables(
        reads, al, bl, seed=12, n_genomes=20, total_bases=2e8, linked=0.0, plasmids=1.0, n_truth_anchors=1000
    )
    ln = s.lengths
    assert {"genome", "replicon"} <= set(ln.columns)
    pl = ln[ln.replicon >= 20]
    assert len(pl) > 0 and pl.genome.max() < 20
    assert (pl.length <= 150_000).all()


def test_chimeras_join_reads_and_keep_features_inside():
    reads, al, bl = _p2_reads(3)
    kw = dict(seed=13, n_genomes=30, total_bases=2e8, linked=0.0, min_frac=0.8, n_truth_anchors=1000)
    a = simulate_tables(reads, al, bl, **kw)
    b = simulate_tables(reads, al, bl, chimera_rate=0.05, **kw)
    assert abs(b.lengths.chimera.mean() - 0.05) < 0.005
    assert b.lengths.length.sum() == a.lengths.length.sum() and len(b.features) == len(a.features)
    m = b.features.merge(b.lengths, on="molecule_id")
    assert (m.start >= 0).all() and (m.end <= m.length).all()
