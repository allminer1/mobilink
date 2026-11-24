"""`mobilink run`: from reads (or contigs) to the net linked test, a per-gene table and a report, in one command.

Steps (each writes its file through a temporary name and is skipped when the file exists, so a stopped run resumes):
  1. lengths                  NAME.lengths.tsv.gz
  2. DIAMOND blastx           NAME.card.tsv.gz, NAME.bacmet.tsv.gz, NAME.mobileog.tsv.gz
  3. features (primary rule)  NAME.features.tsv.gz      ARG = CARD >= 80%, MRG = BacMet2 >= 80%, MGE = mobility core
                                                        of mobileOG-db >= 60%; subject cover >= 80%
  4. host labels (optional)   NAME.kraken.tsv.gz        Kraken2 at confidence 0, all molecules
  5. overlaps (reads only)    NAME.ava.paf.gz           minimap2 all-against-all among the molecules with a call
  6. chimera flags (reads)    NAME.chimeras.tsv
  7. test                     NAME.test.tsv             the net linked test; NAME.test_no_chimeras.tsv as sensitivity
  8. curve, genes, screen     NAME.curve.tsv, NAME.genes.tsv, NAME.fixed_length.tsv (flagged clusters are checked
                              against UniVec; a laboratory plasmid gives NAME.test_no_lab_plasmid.tsv)
  9. report                   NAME.report.html, NAME.run.json

The annotation settings are those of the survey in the paper (DIAMOND long-read mode: -F 15 --range-culling --top 10).
"""

from __future__ import annotations

import gzip
import json
import os
import shutil
import subprocess
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from . import __version__, dbsetup
from . import databases as dbl
from .chimeras import chimera_flags, read_sequences
from .features import RULES, call_features, check_hit_lengths, summarise
from .genes import fixed_length_clusters, gene_table, pick_reads, vector_rule
from .hosts import UNRESOLVED, read_taxonomy, strata_for_molecules
from .io import read_features, read_lengths
from .io.diamond import read_hits
from .loci import read_paf
from .pairs import net_test_pairs
from .stats.nettest import net_w_estimate

DIAMOND_FIELDS = "qseqid sseqid pident length qstart qend sstart send slen qlen evalue bitscore".split()
CURVE_W = (250, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 7500, 10000)


class PipelineError(RuntimeError):
    """A step failed; the message names the step and the log file."""


@dataclass
class RunConfig:
    reads: str
    db_dir: str
    out_dir: str
    name: str | None = None
    threads: int = 8
    kraken_db: str | None = None
    contigs: bool = False
    pairs: tuple = (("ARG", "MGE"), ("MRG", "MGE"), ("ARG", "MRG"))
    w: int = 2000
    units_d: int = 1000
    n_mc: int = 9999
    seed: int = 20261114
    tmp_dir: str | None = None
    keep_hits: bool = True
    tools: dict = field(default_factory=lambda: dict(diamond="diamond", minimap2="minimap2", kraken2="kraken2"))

    def path(self, suffix: str) -> str:
        return os.path.join(self.out_dir, f"{self.name}.{suffix}")


def _exe(tool: str) -> str:
    exe = shutil.which(tool)
    if exe is None:
        raise PipelineError(f"{tool} is not on the PATH (install it, e.g. `conda install -c bioconda "
                            f"{os.path.basename(tool)}`)")
    return exe


def _call(cmd: list, log_path: str, step: str, stdout=None) -> None:
    with open(log_path, "w") as log:
        r = subprocess.run(cmd, stdout=stdout or log, stderr=log, check=False)  # noqa: S603
    if r.returncode != 0:
        with open(log_path) as fh:
            tail = "".join(fh.readlines()[-15:])
        raise PipelineError(f"step '{step}' failed (exit {r.returncode}); log: {log_path}\n{tail}")


def _has_rows(gz_path: str) -> bool:
    with gzip.open(gz_path, "rt") as fh:
        return bool(fh.readline())


def default_name(reads: str) -> str:
    n = os.path.basename(reads)
    for ext in (".gz", ".fastq", ".fq", ".fasta", ".fa", ".fna"):
        if n.endswith(ext):
            n = n[: -len(ext)]
    return n or "sample"


def write_lengths(reads: str, out: str) -> dict:
    """Lengths table of all molecules; returns molecules, bases and N50."""
    n, lens = 0, []
    with gzip.open(out + ".part", "wt") as fo:
        fo.write("molecule_id\tlength\n")
        for rid, seq in read_sequences(reads):
            fo.write(f"{rid}\t{len(seq)}\n")
            lens.append(len(seq))
            n += 1
    if n == 0:
        os.remove(out + ".part")
        raise PipelineError(f"{reads}: no sequence found")
    os.replace(out + ".part", out)
    return length_stats(np.array(lens, dtype=np.int64))


def length_stats(lens: np.ndarray) -> dict:
    s = np.sort(lens)[::-1]
    c = np.cumsum(s)
    return dict(molecules=int(len(s)), bases=int(c[-1]), n50=int(s[np.searchsorted(c, c[-1] / 2)]),
                median_length=int(np.median(s)))


def diamond_blastx(cfg: RunConfig, tag: str, dmnd: str) -> str:
    out = cfg.path(f"{tag}.tsv.gz")
    if os.path.exists(out):
        return out
    tmp = cfg.path(f"{tag}.part.tsv")  # DIAMOND adds .gz
    cmd = [_exe(cfg.tools["diamond"]), "blastx", "-q", cfg.reads, "-d", dmnd, "-F", "15", "--range-culling", "--top",
           "10", "-p", str(cfg.threads), "-f", "6", *DIAMOND_FIELDS, "--compress", "1", "-o", tmp]
    if cfg.tmp_dir:
        cmd += ["--tmpdir", cfg.tmp_dir]
    _call(cmd, cfg.path(f"{tag}.diamond.log"), f"DIAMOND against {tag}")
    os.replace(tmp + ".gz", out)
    return out


def build_features(cfg: RunConfig, db: dict, lengths: pd.DataFrame, hits: dict) -> pd.DataFrame:
    out = cfg.path("features.tsv.gz")
    if os.path.exists(out):
        return read_features(out)
    p = db["paths"]
    classes = dbl.load_bacmet_classes(p["bacmet_classes"])
    lookup = dbl.load_mobileog_lookup(p["mobileog_lookup"])
    core = dbl.load_mobility_core(p["mobility_core"])
    labellers = {"ARG": dbl.label_card, "MRG": lambda d: dbl.label_bacmet(d, classes),
                 "MGE": lambda d: dbl.label_mobileog(d, lookup, core=core)}
    rule = RULES["primary"]
    raw = {}
    for cls, tag in (("ARG", "card"), ("MRG", "bacmet"), ("MGE", "mobileog")):
        raw[cls] = read_hits(hits[tag], rule[cls])
        missing, wrong = check_hit_lengths(raw[cls], lengths)
        if missing or wrong:
            raise PipelineError(f"{cls}: {missing} hit queries missing from the lengths table, {wrong} with a "
                                "different length (are the hit tables from another input file?)")
    f = call_features(raw, rule, labellers)
    f.to_csv(out + ".part", sep="\t", index=False, compression="gzip")
    os.replace(out + ".part", out)
    pd.DataFrame(summarise(f, rule="primary")).to_csv(cfg.path("features.summary.tsv"), sep="\t", index=False)
    return read_features(out)


def kmer_counts(taxid: str, lca: str) -> tuple[int, int]:
    """From the LCA column of a Kraken2 output row (e.g. "562:13 561:4 A:31 0:1"): k-mers that hit any taxon, and
    k-mers that hit the assigned taxon itself. Unclassified (0), ambiguous (A) and the mate separator are skipped."""
    hit = own = 0
    for tok in lca.split():
        t, _, k = tok.partition(":")
        if t not in ("0", "A", "|"):
            hit += int(k)
            if t == taxid:
                own += int(k)
    return hit, own


def kraken_labels(cfg: RunConfig) -> str:
    """Kraken2 at confidence 0 on all molecules -> read_id, status, taxid, hit_kmers, taxon_kmers."""
    out = cfg.path("kraken.tsv.gz")
    if os.path.exists(out):
        return out
    raw = cfg.path("kraken.part.out")
    cmd = [_exe(cfg.tools["kraken2"]), "--db", cfg.kraken_db, "--threads", str(cfg.threads), "--confidence", "0",
           "--report", cfg.path("kraken.report"), "--output", raw, cfg.reads]
    _call(cmd, cfg.path("kraken.log"), "Kraken2")
    with open(raw) as fh, gzip.open(out + ".part", "wt") as fo:
        fo.write("read_id\tstatus\ttaxid\thit_kmers\ttaxon_kmers\n")
        for line in fh:
            status, rid, taxid, _, lca = line.rstrip("\n").split("\t")[:5]
            hit, own = kmer_counts(taxid, lca)
            fo.write(f"{rid}\t{status}\t{taxid}\t{hit}\t{own}\n")
    os.replace(out + ".part", out)
    os.remove(raw)
    return out


def carrier_overlaps(cfg: RunConfig, carriers: set) -> str:
    """minimap2 all-against-all among the carrier molecules; blocks of at least 1 kb, the first 12 PAF columns."""
    out = cfg.path("ava.paf.gz")
    if os.path.exists(out):
        return out
    fa = cfg.path("carriers.part.fa")
    with open(fa, "w") as fo:
        for rid, seq in read_sequences(cfg.reads):
            if rid in carriers:
                fo.write(f">{rid}\n{seq}\n")
    cmd = [_exe(cfg.tools["minimap2"]), "-x", "ava-ont", "-t", str(cfg.threads), fa, fa]
    with open(cfg.path("minimap2.log"), "w") as log, gzip.open(out + ".part", "wt") as fo:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log, text=True)  # noqa: S603
        for line in p.stdout:
            f = line.split("\t", 12)
            if int(f[10]) >= 1000:
                fo.write("\t".join(f[:12]).rstrip("\n") + "\n")
        if p.wait() != 0:
            raise PipelineError(f"step 'minimap2 overlaps' failed (exit {p.returncode}); log: "
                                f"{cfg.path('minimap2.log')}")
    os.replace(out + ".part", out)
    os.remove(fa)
    return out


def vector_hits(cfg: RunConfig, univec: str, picks: pd.DataFrame) -> pd.DataFrame:
    """The picked reads aligned to UniVec (minimap2 map-ont with secondary alignments, as in the survey screen)."""
    fa = cfg.path("vector_reads.part.fa")
    want = set(picks.molecule_id)
    with open(fa, "w") as fo:
        for rid, seq in read_sequences(cfg.reads):
            if rid in want:
                fo.write(f">{rid}\n{seq}\n")
    cmd = [_exe(cfg.tools["minimap2"]), "-x", "map-ont", "-c", "--secondary=yes", "-N", "30", "-t", str(cfg.threads),
           univec, fa]
    with open(cfg.path("vector.log"), "w") as log:
        r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=log, text=True, check=False)  # noqa: S603
    os.remove(fa)
    if r.returncode != 0:
        raise PipelineError(f"step 'UniVec check' failed (exit {r.returncode}); log: {cfg.path('vector.log')}")
    rows = [(f[0], int(f[2]), int(f[3]), f[5], int(f[10])) for f in (ln.split("\t") for ln in r.stdout.splitlines())]
    return pd.DataFrame(rows, columns=["read_id", "q_start", "q_end", "univec_id", "block"])


def curve(lengths, features, pairs, strata, units_d, seed, ws=CURVE_W) -> pd.DataFrame:
    """net_w at a grid of windows (the estimate only): how linkage builds up with distance. 0 is chance at any w."""
    rows = []
    for a, b in pairs:
        fa, fb = features[features["class"] == a], features[features["class"] == b]
        if fa.empty or fb.empty:
            continue
        for w in ws:
            try:
                r = net_w_estimate(lengths, fa, fb, strata=strata, w=w, units_d=units_d, seed=seed)
            except ValueError:
                break
            rows.append(dict(class_a=a, class_b=b, w=w, net_w=r["net_w"], sides_w=r["sides_w"],
                             events_w=r["events_w"]))
    return pd.DataFrame(rows, columns=["class_a", "class_b", "w", "net_w", "sides_w", "events_w"])


def run(cfg: RunConfig, log=print) -> dict:
    """Run all steps; returns the dictionary written to NAME.run.json."""
    from .report import write_report

    cfg.name = cfg.name or default_name(cfg.reads)
    os.makedirs(cfg.out_dir, exist_ok=True)
    db = dbsetup.load(cfg.db_dir)
    for tool in ["diamond"] + ([] if cfg.contigs else ["minimap2"]) + (["kraken2"] if cfg.kraken_db else []):
        _exe(cfg.tools[tool])
    taxonomy = None
    if cfg.kraken_db:
        taxonomy = os.path.join(cfg.kraken_db, "ktaxonomy.tsv")
        if not os.path.exists(taxonomy):
            raise PipelineError(f"{cfg.kraken_db}: no ktaxonomy.tsv (the prebuilt Kraken2 indexes include it)")
    info = dict(mobilink=__version__, started=datetime.now(timezone.utc).isoformat(timespec="seconds"),
                settings={k: v for k, v in asdict(cfg).items() if k != "tools"},
                databases={k: {x: v.get(x) for x in ("version", "source", "fasta_md5", "licence", "cite")}
                           for k, v in db["databases"].items()}, flags=[])

    log("[1/9] lengths")
    lp = cfg.path("lengths.tsv.gz")
    if not os.path.exists(lp):
        write_lengths(cfg.reads, lp)
    lengths = read_lengths(lp)
    info["input"] = length_stats(lengths.length.to_numpy(np.int64))

    hits = {}
    for i, (tag, key) in enumerate((("card", "card"), ("bacmet", "bacmet"), ("mobileog", "mobileog"))):
        log(f"[2/9] DIAMOND blastx against {tag} ({i + 1} of 3)")
        hits[tag] = diamond_blastx(cfg, tag, db["paths"][key])

    log("[3/9] features")
    features = build_features(cfg, db, lengths, hits)
    features = features[features.molecule_id.isin(lengths.molecule_id)]
    info["calls"] = {c: int(n) for c, n in features["class"].value_counts().items()}
    if not cfg.keep_hits:
        for p in hits.values():
            os.remove(p)

    strata = None
    if cfg.kraken_db:
        log("[4/9] host labels (Kraken2)")
        k = pd.read_csv(kraken_labels(cfg), sep="\t", dtype={"read_id": str})
        strata = strata_for_molecules(lengths.molecule_id.to_numpy(), k, read_taxonomy(taxonomy), rank="S",
                                      min_hits=10)
        info["host_unresolved_share"] = float(np.mean(strata == UNRESOLVED))
    else:
        log("[4/9] host labels: skipped (no --kraken-db); the test then assumes one host community")
        info["flags"].append("no host strata: genes that share a host count as linked beyond chance less strictly; "
                             "give --kraken-db for the model of the paper")

    classes = {c for pair in cfg.pairs for c in pair}
    a_classes = {a for a, _ in cfg.pairs}
    paf, exclude = None, None
    if cfg.contigs:
        log("[5/9] overlaps: skipped (contigs; all units of one contig form one clump)")
        log("[6/9] chimera flags: skipped (contigs)")
    else:
        log("[5/9] overlaps among molecules with a call (minimap2)")
        carriers = set(features.molecule_id[features["class"].isin(classes)])
        pp = carrier_overlaps(cfg, carriers)
        log("[6/9] chimera flags")
        cp = cfg.path("chimeras.tsv")
        if not os.path.exists(cp):
            if _has_rows(pp):
                chim = chimera_flags(paf=pp, tmp_dir=cfg.tmp_dir)
            else:  # no overlap at all (shallow sample)
                chim = pd.DataFrame(columns=["read_id", "overlap_assessable", "n_breakpoints", "flag_overlap",
                                             "n_adapters", "flag_adapter", "flag"])
            chim.to_csv(cp, sep="\t", index=False)
        chim = pd.read_csv(cp, sep="\t", dtype={"read_id": str})
        exclude = set(chim.read_id[chim.flag.astype(str).str.lower().isin(["true", "1"])])
        info["chimeras"] = dict(judged=int(len(chim)), flagged=int(len(exclude)))
        paf = read_paf(pp, min_block=1000, keep=set(features.molecule_id[features["class"].isin(a_classes)]))

    log("[7/9] net linked test")
    kw = dict(strata=strata, paf=paf, loci="molecule" if cfg.contigs else None, w=cfg.w, units_d=cfg.units_d,
              seed=cfg.seed, n_mc=cfg.n_mc)
    tp = cfg.path("test.tsv")
    if not os.path.exists(tp):
        net_test_pairs(lengths, features, list(cfg.pairs), **kw).to_csv(tp, sep="\t", index=False)
    test = pd.read_csv(tp, sep="\t")
    sens = None
    if exclude:
        sp = cfg.path("test_no_chimeras.tsv")
        if not os.path.exists(sp):
            net_test_pairs(lengths, features, list(cfg.pairs), exclude=exclude, oe=False, **kw).to_csv(
                sp, sep="\t", index=False)
        sens = pd.read_csv(sp, sep="\t")

    log("[8/9] distance curve, per-gene table, fixed-length screen")
    cvp = cfg.path("curve.tsv")
    if not os.path.exists(cvp):
        curve(lengths, features, cfg.pairs, strata, cfg.units_d, cfg.seed).to_csv(cvp, sep="\t", index=False)
    cv = pd.read_csv(cvp, sep="\t")
    gp = cfg.path("genes.tsv")
    if not os.path.exists(gp):
        ann = None
        if "aro_index" in db["paths"]:
            ai = pd.read_csv(db["paths"]["aro_index"], sep="\t", dtype=str)
            ann = ai.rename(columns={"ARO Accession": "label", "AMR Gene Family": "gene_family",
                                     "Drug Class": "drug_class", "Resistance Mechanism": "mechanism"})[
                ["label", "gene_family", "drug_class", "mechanism"]].drop_duplicates("label")
        parts = [gene_table(lengths, features, a, b, strata=strata, w=cfg.w, units_d=cfg.units_d, seed=cfg.seed,
                            annotation=ann if a == "ARG" else None) for a, b in cfg.pairs]
        parts = [p for p in parts if not p.empty]
        (pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=["class_a", "class_b", "gene"])
         ).to_csv(gp, sep="\t", index=False)
    genes = pd.read_csv(gp, sep="\t")
    fixed, clean = pd.DataFrame(), None
    if not cfg.contigs:
        fp = cfg.path("fixed_length.tsv")
        if not os.path.exists(fp):
            fixed = fixed_length_clusters(lengths, features, "ARG")
            univec = db["paths"].get("univec")
            if univec and fixed.flag.any():
                picks = pick_reads(lengths, features, fixed)
                fixed = fixed.merge(vector_rule(picks, vector_hits(cfg, univec, picks)), on="gene", how="left")
            fixed.to_csv(fp, sep="\t", index=False)
        fixed = pd.read_csv(fp, sep="\t")
        lab = set()
        for r in fixed[fixed.flag].itertuples():
            head = (f"{r.gene}: {r.carriers} carrier reads, {r.share_near_median:.0%} of them within 3% of "
                    f"{r.median_length / 1000:.2f} kb, the pattern of a small plasmid sequenced whole. ")
            if "laboratory_plasmid" not in fixed.columns:
                info["flags"].append(head + "UniVec is not installed, so the reads were not checked; align a few of "
                                     "them to a vector database before you interpret ARG results.")
            elif str(r.laboratory_plasmid).lower() == "true":
                lab.add(r.gene)
                info["flags"].append(
                    head + f"The reads match vector sequences ({r.median_covered:.0%} of the read, "
                    f"{r.median_sources:g} UniVec records at the median): a laboratory plasmid. The table "
                    "'without laboratory plasmids' leaves out all reads that carry this gene; use it.")
            else:
                info["flags"].append(head + f"The reads do not look like a cloning vector ({r.median_covered:.0%} "
                                     "of the read matches UniVec): probably a natural small plasmid; kept.")
        if lab:
            drop = set(features.molecule_id[(features["class"] == "ARG") & features.gene.isin(lab)])
            info["laboratory_plasmid"] = dict(genes=sorted(lab), reads_removed=len(drop))
            vp = cfg.path("test_no_lab_plasmid.tsv")
            if not os.path.exists(vp):
                net_test_pairs(lengths, features, list(cfg.pairs), exclude=drop, oe=False, **kw).to_csv(
                    vp, sep="\t", index=False)
            clean = pd.read_csv(vp, sep="\t")

    info["finished"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    info["files"] = sorted(f for f in os.listdir(cfg.out_dir) if f.startswith(cfg.name + ".") and ".part" not in f)
    log("[9/9] report")
    write_report(cfg.path("report.html"), cfg.name, info, test, sens, cv, genes, fixed, w=cfg.w, contigs=cfg.contigs,
                 clean=clean)
    with open(cfg.path("run.json"), "w") as fh:
        json.dump(info, fh, indent=2, default=str)
    return info
