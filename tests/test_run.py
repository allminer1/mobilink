"""`mobilink db` and `mobilink run`: helpers without external tools, and the whole chain when DIAMOND and minimap2 are
installed (skipped otherwise)."""

import gzip
import io
import json
import shutil
import tarfile

import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from mobilink import dbsetup
from mobilink.cli import main
from mobilink.genes import fixed_length_clusters, gene_table
from mobilink.pipeline import default_name, kmer_counts, length_stats
from mobilink.report import verdict, write_report

AA = "ACDEFGHIKLMNPQRSTVWY"
CODON = dict(zip(AA, "GCT TGT GAT GAA TTT GGT CAT ATT AAA CTG ATG AAT CCG CAG CGT TCT ACC GTG TGG TAT".split(),
                 strict=True))
COMP = str.maketrans("ACGT", "TGCA")


def test_packaged_tables_are_installed():
    assert pd.read_csv(dbsetup.data_file("bacmet2_exp_classes.tsv"), sep="\t").shape[0] == 753
    assert pd.read_csv(dbsetup.data_file("mobility_core.tsv.gz"), sep="\t", usecols=["sseqid"]).shape[0] == 109_325


def test_kmer_counts_skip_unclassified_and_ambiguous():
    assert kmer_counts("562", "562:13 561:4 A:31 0:1 |:| 562:2") == (19, 15)
    assert kmer_counts("0", "0:40") == (0, 0)


def test_length_stats_and_default_name():
    s = length_stats(np.array([1000, 2000, 3000, 4000]))
    assert (s["molecules"], s["bases"], s["n50"]) == (4, 10_000, 3000)
    assert default_name("/x/sample1.fastq.gz") == "sample1"


def test_mobileog_lookup_restores_fields_cut_by_diamond(tmp_path):
    fa = tmp_path / "m.faa"
    fa.write_text(">mobileOG_1|int|P1|integration/excision|NA|Plasmid RefSeq|Keyword Search\nMK\n"
                  ">mobileOG_2|tra|P2|transfer|NA|Multiple|Manual\nMK\n")
    out = tmp_path / "lookup.tsv.gz"
    assert dbsetup.mobileog_lookup(str(fa), str(out)) == 2
    t = pd.read_csv(out, sep="\t")
    assert t.sseqid[0] == "mobileOG_1|int|P1|integration/excision|NA|Plasmid"
    assert list(t.category) == ["integration/excision", "transfer"] and t.evidence[0] == "Keyword Search"


def test_card_archive_gives_fasta_index_and_version(tmp_path):
    arc = tmp_path / "card.tar.bz2"
    with tarfile.open(arc, "w:bz2") as tf:
        for name, text in ((dbsetup.CARD_MEMBER, ">gb|X|ARO:1|geneA [sp]\nMK\n"), ("aro_index.tsv", "ARO Accession\n"),
                           ("card.json", json.dumps({"_version": "9.9.9"}))):
            data = text.encode()
            ti = tarfile.TarInfo("./" + name)
            ti.size = len(data)
            tf.addfile(ti, io.BytesIO(data))
    fasta, index, version = dbsetup.card_from_archive(str(arc), str(tmp_path))
    assert version == "9.9.9" and index.endswith("aro_index.tsv") and dbsetup.count_fasta(fasta) == 1


def test_load_refuses_a_directory_without_manifest(tmp_path):
    with pytest.raises(dbsetup.DbError, match="not a mobilink database directory"):
        dbsetup.load(str(tmp_path))


def _features(rows):
    return pd.DataFrame(rows, columns=["molecule_id", "start", "end", "class", "gene", "label"]).assign(scov=100.0)


def test_gene_table_counts_whole_neighbours_within_the_window():
    ln = pd.DataFrame({"molecule_id": ["m1", "m2", "m3"], "length": [10_000, 10_000, 10_000]})
    f = _features([("m1", 3000, 3900, "ARG", "argA", "ARO:1"),
                   ("m1", 4100, 5000, "MGE", "intI1", "integration/excision"),
                   ("m2", 3000, 3900, "ARG", "argA", "ARO:1"), ("m2", 6500, 7400, "MGE", "traA", "transfer"),
                   ("m3", 3000, 3900, "ARG", "argB", "ARO:2")])
    t = gene_table(ln, f, "ARG", "MGE", w=2000).set_index("gene")
    assert t.loc["argA", "calls"] == 2 and t.loc["argA", "calls_with_b"] == 1  # the far end on m2 is 3.5 kb away
    assert t.loc["argA", "top_b_label"] == "integration/excision" and t.loc["argA", "top_b_gene"] == "intI1"
    assert t.loc["argB", "calls_with_b"] == 0 and np.isnan(t.loc["argA", "net_w"])  # fewer than 30 calls


def test_fixed_length_screen_flags_one_read_length_only():
    g = np.random.default_rng(1)
    plasmid = [(f"p{i}", int(6070 + g.integers(-60, 60))) for i in range(40)]
    spread = [(f"s{i}", int(g.integers(3000, 30_000))) for i in range(40)]
    few = [(f"f{i}", 6000) for i in range(10)]
    ln = pd.DataFrame(plasmid + spread + few, columns=["molecule_id", "length"])
    f = _features([(m, 100, 900, "ARG", "aph", "ARO:1") for m, _ in plasmid]
                  + [(m, 100, 900, "ARG", "tetW", "ARO:2") for m, _ in spread]
                  + [(m, 100, 900, "ARG", "rare", "ARO:3") for m, _ in few])
    t = fixed_length_clusters(ln, f).set_index("gene")
    assert bool(t.loc["aph", "flag"]) and not bool(t.loc["tetW", "flag"]) and "rare" not in t.index


def test_verdict_wording_follows_estimate_and_p_value():
    base = dict(class_a="ARG", class_b="MGE", net_w_lo=0.05, net_w_hi=0.15)
    assert "are linked: 10.0% of the flanks" in verdict(dict(base, net_w=0.10, net_w_p=0.001), 2000)
    assert "p < 0.001" in verdict(dict(base, net_w=0.10, net_w_p=0.0001), 2000)
    assert "includes 0 here" in verdict(dict(base, net_w=0.03, net_w_lo=-0.01, net_w_p=0.001), 2000)
    assert "no linkage shown" in verdict(dict(base, net_w=0.01, net_w_p=0.6), 2000)
    assert "no estimate" in verdict(dict(base, net_w=np.nan, net_w_p=np.nan, flags="class empty"), 2000)


def test_report_is_one_self_contained_file(tmp_path):
    test = pd.DataFrame([dict(class_a="ARG", class_b="MGE", net_w=0.1, net_w_lo=0.05, net_w_hi=0.15, net_w_p=0.001,
                              net_w_q=0.002, units=120, sides_w=150, net_w_O=20, net_w_E=6.5, net_mol=0.2, flags="")])
    cv = pd.DataFrame(dict(class_a="ARG", class_b="MGE", w=[500, 1000, 2000, 5000], net_w=[0.02, 0.06, 0.1, 0.12],
                           sides_w=[200, 180, 150, 90], events_w=[5, 12, 20, 22]))  # the last window is not drawn
    info = dict(mobilink="x", finished="2026-11-24T00:00:00", input=dict(molecules=10, bases=5e4, n50=6000),
                calls=dict(ARG=3), flags=["a <note>"], databases={}, files=["s.test.tsv"])
    out = tmp_path / "r.html"
    write_report(str(out), "s<1>", info, test, None, cv, pd.DataFrame(), pd.DataFrame())
    h = out.read_text()
    assert h.count("<circle") == 3
    assert "<svg" in h and "s&lt;1&gt;" in h and "a &lt;note&gt;" in h and "<script" not in h and "http" not in h


def _protein(g, n=300):
    return "M" + "".join(g.choice(list(AA), n - 1))


def _dna(p):
    return "".join(CODON[a] for a in p) + "TAA"


@pytest.mark.skipif(shutil.which("diamond") is None or shutil.which("minimap2") is None,
                    reason="needs DIAMOND and minimap2 on the PATH")
def test_db_and_run_end_to_end(tmp_path):
    g = np.random.default_rng(7)
    core = pd.read_csv(dbsetup.data_file("mobility_core.tsv.gz"), sep="\t", usecols=["sseqid"]).sseqid[:2].tolist()
    arg = {"gb|T1|ARO:3000001|argA [test]": _protein(g), "gb|T2|ARO:3000002|argB [test]": _protein(g)}
    mrg = {"BAC0003|acn|tr|O53166|TEST": _protein(g), "BAC0004|acr3|tr|B5LX01|TEST": _protein(g)}
    mge = {core[0]: _protein(g), core[1]: _protein(g),
           "mobileOG_999999999|foo|X1|phage|NA|Multiple|Manual": _protein(g)}
    for name, d in (("card", arg), ("bacmet", mrg), ("mog", mge)):
        (tmp_path / f"{name}.faa").write_text("".join(f">{h}\n{s}\n" for h, s in d.items()))
    a, m, e = list(arg.values()), list(mrg.values()), list(mge.values())
    genomes = []
    # every ARG has a core MGE gene 200 bp from one flank and none near the other, so the truth is close to 0.5;
    # the MRG and the non-core gene sit far from everything
    for i in range(6):
        rnd = ["".join(g.choice(list("ACGT"), n)) for n in (15_000, 200, 12_000, 200, 12_000, 9_000, 9_000)]
        genomes.append(rnd[0] + _dna(a[i % 2]) + rnd[1] + _dna(e[i % 2]) + rnd[2] + _dna(a[(i + 1) % 2]) + rnd[3]
                       + _dna(e[(i + 1) % 2]) + rnd[4] + _dna(m[i % 2]) + rnd[5] + _dna(e[2]) + rnd[6])
    reads = tmp_path / "reads.fa.gz"
    with gzip.open(reads, "wt") as fo:
        for k in range(3000):
            s = genomes[k % 6]
            ln = int(g.integers(3000, 9000))
            st = int(g.integers(0, len(s) - ln))
            r = s[st:st + ln]
            fo.write(f">r{k}\n{r if k % 2 else r.translate(COMP)[::-1]}\n")
    run = CliRunner()
    r = run.invoke(main, ["db", "--out", str(tmp_path / "db"), "--card", str(tmp_path / "card.faa"), "--bacmet",
                          str(tmp_path / "bacmet.faa"), "--mobileog", str(tmp_path / "mog.faa"), "--threads", "2"])
    assert r.exit_code == 0, r.output
    args = ["run", str(reads), "--db", str(tmp_path / "db"), "--out", str(tmp_path / "out"), "--name", "toy",
            "--threads", "2", "--n-mc", "199"]
    r = run.invoke(main, args)
    assert r.exit_code == 0, r.output
    out = tmp_path / "out"
    f = pd.read_csv(out / "toy.features.tsv.gz", sep="\t")
    assert {"ARG", "MRG", "MGE", "MGE_other"} <= set(f["class"])
    t = pd.read_csv(out / "toy.test.tsv", sep="\t").set_index(["class_a", "class_b"])
    assert 0.4 < t.loc[("ARG", "MGE"), "net_w"] < 0.6 and t.loc[("ARG", "MGE"), "net_w_p"] < 0.05
    assert not t.loc[("MRG", "MGE"), "net_w"] > 0.2
    genes = pd.read_csv(out / "toy.genes.tsv", sep="\t")
    top = genes[(genes.class_a == "ARG") & (genes.class_b == "MGE")].set_index("gene")
    assert top.loc["argA", "share_with_b"] > 0.5
    assert "<svg" in (out / "toy.report.html").read_text()
    assert json.loads((out / "toy.run.json").read_text())["calls"]["ARG"] > 100
    r = run.invoke(main, args)  # a second call finds every file and only rewrites the report
    assert r.exit_code == 0, r.output
    r = run.invoke(main, ["run", str(tmp_path / "card.faa"), "--db", str(tmp_path / "nodb"), "--out", str(out)])
    assert r.exit_code != 0
