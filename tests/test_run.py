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
from mobilink.genes import fixed_length_clusters, gene_table, pick_reads, vector_rule
from mobilink.pipeline import default_name, kmer_counts, length_stats
from mobilink.power import depth_projection
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


def test_vector_rule_separates_a_cloning_vector_from_a_natural_plasmid():
    ln = pd.DataFrame({"molecule_id": [f"v{i}" for i in range(40)] + [f"n{i}" for i in range(40)], "length": 6000})
    f = _features([(f"v{i}", 100, 900, "ARG", "aph", "ARO:1") for i in range(40)]
                  + [(f"n{i}", 100, 900, "ARG", "qnrB", "ARO:2") for i in range(40)])
    picks = pick_reads(ln, f, fixed_length_clusters(ln, f))
    assert picks.groupby("gene").size().to_dict() == {"aph": 4, "qnrB": 4}
    hits = [(m, 1000, 2600, "gnl|uv|J01749.1:1-4361", 1600) for m in picks.molecule_id[picks.gene == "aph"]]
    hits += [(m, 3000, 4500, "gnl|uv|X02513.1:1-1500", 1500) for m in picks.molecule_id[picks.gene == "aph"]]
    hits += [(m, 200, 1100, "gnl|uv|J01749.1:1-4361", 900) for m in picks.molecule_id[picks.gene == "qnrB"]]
    t = vector_rule(picks, pd.DataFrame(hits, columns=["read_id", "q_start", "q_end", "univec_id", "block"]))
    t = t.set_index("gene")
    assert bool(t.loc["aph", "laboratory_plasmid"]) and not bool(t.loc["qnrB", "laboratory_plasmid"])
    assert abs(t.loc["aph", "median_covered"] - 3100 / 6000) < 1e-9 and t.loc["qnrB", "median_sources"] == 1


def test_verdict_wording_follows_estimate_and_p_value():
    base = dict(class_a="ARG", class_b="MGE", net_w_lo=0.05, net_w_hi=0.15)
    assert "are linked: 10.0% of the flanks" in verdict(dict(base, net_w=0.10, net_w_p=0.001), 2000)
    assert "p < 0.001" in verdict(dict(base, net_w=0.10, net_w_p=0.0001), 2000)
    assert "includes 0 here" in verdict(dict(base, net_w=0.03, net_w_lo=-0.01, net_w_p=0.001), 2000)
    assert "no linkage shown" in verdict(dict(base, net_w=0.01, net_w_p=0.6), 2000)
    assert "no estimate" in verdict(dict(base, net_w=np.nan, net_w_p=np.nan, flags="class empty"), 2000)


def test_depth_projection_uses_the_test_for_depth_and_the_standard_error_for_the_interval():
    t = pd.DataFrame([dict(class_a="ARG", class_b="MGE", net_w=0.1, net_w_se=0.02, net_w_z=20.0),
                      dict(class_a="MRG", class_b="MGE", net_w=np.nan, net_w_se=np.nan, net_w_z=np.nan)])
    p = depth_projection(t, bases=3e9).set_index("class_a")
    assert p.loc["ARG", "halfwidth"] == pytest.approx(1.96 * 0.02, rel=1e-3)
    assert p.loc["ARG", "test_sd"] == pytest.approx(0.005)
    assert p.loc["ARG", "detectable_net"] == pytest.approx((1.96 + 0.8416) * 0.005, rel=1e-3)
    assert p.loc["ARG", "gbp_to_detect_0.01"] == pytest.approx(3 * ((1.96 + 0.8416) * 0.005 / 0.01) ** 2, rel=1e-3)
    assert p.loc["ARG", "gbp_to_detect_0.01"] == pytest.approx(4 * p.loc["ARG", "gbp_to_detect_0.02"])
    assert np.isnan(p.loc["MRG", "halfwidth"]) and np.isnan(p.loc["MRG", "gbp_to_detect_0.05"])
    assert "gbp" not in depth_projection(t).columns  # contigs: no depth columns
    no_estimate = pd.DataFrame([dict(class_a="ARG", class_b="MRG", net_w=np.nan, net_w_se=0.0003, net_w_z=np.nan)])
    assert np.isnan(depth_projection(no_estimate, bases=1e9).halfwidth.iloc[0])  # no half-width without an estimate
    wide = pd.DataFrame([dict(class_a="ARG", class_b="MGE", net_w=0.05, net_w_se=0.014, net_w_z=3.4, net_w_lo=-0.04,
                              net_w_hi=0.14)])
    assert depth_projection(wide).halfwidth.iloc[0] == pytest.approx(0.09)  # the reported interval, not 1.96 * se


def test_every_alignment_is_kept_unless_the_calling_of_the_paper_is_asked_for(tmp_path, monkeypatch):
    from mobilink import pipeline

    seen = []

    def fake_call(cmd, log_path, step, stdout=None):
        seen.append(cmd)
        with gzip.open(cmd[cmd.index("-o") + 1] + ".gz", "wt"):
            pass

    monkeypatch.setattr(pipeline, "_call", fake_call)
    monkeypatch.setattr(pipeline, "_exe", lambda tool: tool)
    for name, one in (("a", False), ("b", True)):
        cfg = pipeline.RunConfig(reads="r.fq", db_dir="db", out_dir=str(tmp_path), name=name, one_alignment=one)
        pipeline.diamond_blastx(cfg, "card", "card.dmnd")
    assert seen[0][seen[0].index("--max-hsps") + 1] == "0" and "--max-hsps" not in seen[1]


def test_report_is_one_self_contained_file(tmp_path):
    test = pd.DataFrame([dict(class_a="ARG", class_b="MGE", net_w=0.1, net_w_lo=0.05, net_w_hi=0.15, net_w_p=0.001,
                              net_w_q=0.002, units=120, sides_w=150, net_w_O=20, net_w_E=6.5, net_mol=0.2, flags="")])
    cv = pd.DataFrame(dict(class_a="ARG", class_b="MGE", w=[500, 1000, 2000, 5000], net_w=[0.02, 0.06, 0.1, 0.12],
                           sides_w=[200, 180, 150, 90], events_w=[5, 12, 20, 22]))  # the last window is not drawn
    info = dict(mobilink="x", finished="2026-11-24T00:00:00", input=dict(molecules=10, bases=5e4, n50=6000),
                calls=dict(ARG=3), flags=["a <note>"], databases={}, files=["s.test.tsv"])
    out = tmp_path / "r.html"
    test["net_w_se"], test["net_w_z"] = 0.025, 10.0
    pw = depth_projection(test, bases=2e9)
    write_report(str(out), "s<1>", info, test, None, cv, pd.DataFrame(), pd.DataFrame(), power=pw)
    h = out.read_text()
    assert "What this depth can show" in h and "±0.050" in h and "0.028" in h and "reached" in h and "16 Gbp" in h
    assert h.count("<circle") == 3
    assert "<svg" in h and "s&lt;1&gt;" in h and "a &lt;note&gt;" in h and "<script" not in h and "http" not in h


def _protein(g, n=300):
    return "M" + "".join(g.choice(list(AA), n - 1))


def _dna(p):
    return "".join(CODON[a] for a in p) + "TAA"


@pytest.mark.skipif(shutil.which("diamond") is None or shutil.which("minimap2") is None,
                    reason="needs DIAMOND and minimap2 on the PATH")
def test_fetch_resumes_a_download_that_is_cut(tmp_path):
    """A server that closes the connection halfway: the download is resumed and the whole file arrives."""
    import threading
    from http.server import BaseHTTPRequestHandler, HTTPServer

    body, seen = bytes(range(256)) * 400, []

    class H(BaseHTTPRequestHandler):
        def do_GET(self):  # noqa: N802
            rng = self.headers.get("Range")
            seen.append(rng)
            start = int(rng.split("=")[1].rstrip("-")) if rng else 0
            self.send_response(206 if rng else 200)
            self.send_header("Content-Length", str(len(body) - start))
            self.end_headers()
            self.wfile.write(body[start:] if rng else body[: len(body) // 2])  # the first answer is cut

        def log_message(self, *a):
            pass

    srv = HTTPServer(("127.0.0.1", 0), H)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    try:
        url, dest, log = f"http://127.0.0.1:{srv.server_port}/f", tmp_path / "f", []
        dbsetup.fetch(url, str(dest), log.append, pause=0)
        assert dest.read_bytes() == body and seen == [None, f"bytes={len(body) // 2}-"]
        assert any("trying again" in x for x in log) and not (tmp_path / "f.part").exists()
        body_full, body = body, body[:10]  # a server that never sends the whole file: an error, and no file is left
        seen.clear()

        class Cut(H):
            def do_GET(self):  # noqa: N802
                self.send_response(200)
                self.send_header("Content-Length", str(len(body_full)))
                self.end_headers()
                self.wfile.write(body)

        srv.RequestHandlerClass = Cut
        with pytest.raises(dbsetup.DbError, match="download failed"):
            dbsetup.fetch(url, str(tmp_path / "g"), lambda m: None, tries=2, pause=0)
        assert not (tmp_path / "g").exists() and not (tmp_path / "g.part").exists()
    finally:
        srv.shutdown()


def test_db_and_run_end_to_end(tmp_path):
    g = np.random.default_rng(7)
    core = pd.read_csv(dbsetup.data_file("mobility_core.tsv.gz"), sep="\t", usecols=["sseqid"]).sseqid[:2].tolist()
    arg = {"gb|T1|ARO:3000001|argA [test]": _protein(g), "gb|T2|ARO:3000002|argB [test]": _protein(g),
           "gb|T3|ARO:3000003|argC [test]": _protein(g)}
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
        # a laboratory plasmid: argC next to two vector segments, sequenced whole 60 times (the gene sits more than
        # 1 kb from the read ends; nearer units are left out of the estimate by the margin rule)
        vec = ["".join(g.choice(list("ACGT"), 1500)) for _ in range(2)]
        plasmid = "".join(g.choice(list("ACGT"), 1300)) + _dna(a[2]) + vec[0] + vec[1] + "".join(
            g.choice(list("ACGT"), 1100))
        for k in range(60):
            r = plasmid[int(g.integers(0, 30)):]
            fo.write(f">p{k}\n{r if k % 2 else r.translate(COMP)[::-1]}\n")
    (tmp_path / "uv.fa").write_text(f">gnl|uv|V1.1:1-1500 test vector one\n{vec[0]}\n"
                                    f">gnl|uv|V2.1:1-1500 test vector two\n{vec[1]}\n")
    run = CliRunner()
    r = run.invoke(main, ["db", "--out", str(tmp_path / "db"), "--card", str(tmp_path / "card.faa"), "--bacmet",
                          str(tmp_path / "bacmet.faa"), "--mobileog", str(tmp_path / "mog.faa"), "--univec",
                          str(tmp_path / "uv.fa"), "--threads", "2"])
    assert r.exit_code == 0, r.output
    args = ["run", str(reads), "--db", str(tmp_path / "db"), "--out", str(tmp_path / "out"), "--name", "toy",
            "--threads", "2", "--n-mc", "199"]
    r = run.invoke(main, args)
    assert r.exit_code == 0, r.output
    out = tmp_path / "out"
    f = pd.read_csv(out / "toy.features.tsv.gz", sep="\t")
    assert {"ARG", "MRG", "MGE", "MGE_other"} <= set(f["class"])
    t = pd.read_csv(out / "toy.test.tsv", sep="\t").set_index(["class_a", "class_b"])
    assert 0.25 < t.loc[("ARG", "MGE"), "net_w"] < 0.6 and t.loc[("ARG", "MGE"), "net_w_p"] < 0.05
    fx = pd.read_csv(out / "toy.fixed_length.tsv", sep="\t").set_index("gene")
    assert bool(fx.loc["argC", "flag"]) and bool(fx.loc["argC", "laboratory_plasmid"])
    assert not fx.drop(index="argC").flag.any()
    c = pd.read_csv(out / "toy.test_no_lab_plasmid.tsv", sep="\t").set_index(["class_a", "class_b"])
    assert c.loc[("ARG", "MGE"), "excluded_molecules"] == 60
    assert 0.4 < c.loc[("ARG", "MGE"), "net_w"] < 0.6
    assert c.loc[("ARG", "MGE"), "net_w"] > t.loc[("ARG", "MGE"), "net_w"] + 0.02  # the plasmid diluted the estimate
    assert "Without laboratory plasmids" in (out / "toy.report.html").read_text()
    assert not t.loc[("MRG", "MGE"), "net_w"] > 0.2
    genes = pd.read_csv(out / "toy.genes.tsv", sep="\t")
    top = genes[(genes.class_a == "ARG") & (genes.class_b == "MGE")].set_index("gene")
    assert top.loc["argA", "share_with_b"] > 0.5
    assert "<svg" in (out / "toy.report.html").read_text()
    pw = pd.read_csv(out / "toy.power.tsv", sep="\t")
    assert len(pw) == len(t) and pw["gbp_to_detect_0.05"].notna().any()
    assert json.loads((out / "toy.run.json").read_text())["calls"]["ARG"] > 100
    r = run.invoke(main, args)  # a second call finds every file and only rewrites the report
    assert r.exit_code == 0, r.output
    r = run.invoke(main, ["run", str(tmp_path / "card.faa"), "--db", str(tmp_path / "nodb"), "--out", str(out)])
    assert r.exit_code != 0

    # contigs are searched in windows: a gene that occurs four times on one contig is called four times, once each
    # (also the copy in the overlap of two windows), at the right places
    pos, parts, at = [8_000, 27_000, 66_000, 101_000], [], 0
    for q in pos:
        parts.append("".join(g.choice(list("ACGT"), q - at)) + _dna(e[0]))
        at = q + len(_dna(e[0]))
    big = "".join(parts) + _dna(a[0]).join("".join(g.choice(list("ACGT"), 9_000)) for _ in range(2))
    small = "".join(g.choice(list("ACGT"), 2_000)) + _dna(e[1]) + "".join(g.choice(list("ACGT"), 2_000))
    (tmp_path / "contigs.fa").write_text(f">big one\n{big}\n>small\n{small}\n")
    r = run.invoke(main, ["run", str(tmp_path / "contigs.fa"), "--contigs", "--db", str(tmp_path / "db"), "--out",
                          str(tmp_path / "outc"), "--name", "c", "--threads", "2", "--n-mc", "199"])
    assert r.exit_code == 0, r.output
    fc = pd.read_csv(tmp_path / "outc" / "c.features.tsv.gz", sep="\t")
    m4 = fc[(fc.molecule_id == "big") & (fc["class"] == "MGE")].sort_values("start")
    assert len(m4) == 4 and all(abs(x - y) <= 3 for x, y in zip(m4.start, pos, strict=True))
    assert len(fc[(fc.molecule_id == "big") & (fc["class"] == "ARG")]) == 1
    assert len(fc[(fc.molecule_id == "small") & (fc["class"] == "MGE")]) == 1
    assert not (tmp_path / "outc" / "c.windows.part.fa").exists()
