"""Command-line entry point. Subcommands are added as their modules are implemented (P1 steps 1.9-1.16)."""

import click
import pandas as pd

from . import __version__


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__, prog_name="linkstat")
def main():
    """Chance-corrected co-localization of gene classes on long reads and contigs."""


@main.command()
@click.option(
    "--lengths",
    "lengths_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Lengths table (molecule_id or read_id, length).",
)
@click.option(
    "--card", type=click.Path(exists=True, dir_okay=False), help="DIAMOND hits against CARD protein homologs."
)
@click.option("--bacmet", type=click.Path(exists=True, dir_okay=False), help="DIAMOND hits against BacMet2 EXP.")
@click.option(
    "--bacmet-classes",
    type=click.Path(exists=True, dir_okay=False),
    help="BacMet class table (BacMet_ID, Accession, Gene_name, metals, is_mrg).",
)
@click.option("--mobileog", type=click.Path(exists=True, dir_okay=False), help="DIAMOND hits against mobileOG-db.")
@click.option(
    "--mobileog-lookup",
    type=click.Path(exists=True, dir_okay=False),
    help="mobileOG lookup (sseqid as cut by DIAMOND, category, ...).",
)
@click.option(
    "--rule",
    "rules",
    multiple=True,
    default=("primary", "strict", "relaxed", "mgeperm"),
    show_default=True,
    type=click.Choice(["primary", "strict", "relaxed", "mgeperm"]),
    help="Threshold rule(s).",
)
@click.option(
    "--mobility-core",
    type=click.Path(exists=True, dir_okay=False),
    help="Mobility-core subject list (column sseqid). If given, class MGE = calls whose subject is in the core "
    "(others: MGE_other); otherwise MGE = the four D3 categories (others: MGE_std).",
)
@click.option("--out-prefix", required=True, help="Writes PREFIX.<rule>.features.tsv.gz and PREFIX.summary.tsv.")
def prep(lengths_path, card, bacmet, bacmet_classes, mobileog, mobileog_lookup, rules, mobility_core, out_prefix):
    """Build feature tables (one call per ORF and database) from DIAMOND hits."""
    from . import databases as db
    from .features import RULES, call_features, check_hit_lengths, summarise
    from .io import read_lengths
    from .io.diamond import read_hits

    if bacmet and not bacmet_classes:
        raise click.UsageError("--bacmet needs --bacmet-classes")
    if mobileog and not mobileog_lookup:
        raise click.UsageError("--mobileog needs --mobileog-lookup")
    lengths = read_lengths(lengths_path, drop_empty=False)
    src = {"ARG": card, "MRG": bacmet, "MGE": mobileog}
    labellers = {"ARG": db.label_card}
    if bacmet:
        classes = db.load_bacmet_classes(bacmet_classes)
        labellers["MRG"] = lambda d: db.label_bacmet(d, classes)
    if mobility_core and not mobileog:
        raise click.UsageError("--mobility-core needs --mobileog")
    if mobileog:
        lookup = db.load_mobileog_lookup(mobileog_lookup)
        core = db.load_mobility_core(mobility_core) if mobility_core else None
        labellers["MGE"] = lambda d: db.label_mobileog(d, lookup, core=core)
    raw = {}
    for cls, path in src.items():
        if path is None:
            continue
        raw[cls] = read_hits(path, min(RULES[r][cls] for r in rules))
        missing, wrong = check_hit_lengths(raw[cls], lengths)
        if missing or wrong:
            raise click.ClickException(
                f"{cls}: {missing} hit queries missing from the lengths table, {wrong} with a different length"
            )
    summary = []
    for rule in rules:
        f = call_features(raw, RULES[rule], labellers)
        f.to_csv(f"{out_prefix}.{rule}.features.tsv.gz", sep="\t", index=False, compression="gzip")
        summary += summarise(f, rule=rule)
    s = pd.DataFrame(summary)
    s.to_csv(f"{out_prefix}.summary.tsv", sep="\t", index=False)
    click.echo(s.to_string(index=False))


@main.command()
@click.option(
    "--read-lengths",
    "read_lengths_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Lengths table whose `length` column is the read-length distribution to draw from.",
)
@click.option(
    "--a-lengths",
    type=click.Path(exists=True, dir_okay=False),
    help="One length per line for class A genes (default: log-normal, median 900 bp).",
)
@click.option(
    "--b-lengths",
    type=click.Path(exists=True, dir_okay=False),
    help="One length per line for class B genes (default: log-normal, median 1,000 bp).",
)
@click.option("--a-class", default="ARG", show_default=True)
@click.option("--b-class", default="MGE", show_default=True)
@click.option("--n-genomes", default=60, show_default=True, type=int)
@click.option("--total-bases", default=1.5e9, show_default=True, type=float)
@click.option("--linked", default=0.3, show_default=True, type=float, help="Share of A genes placed next to a B gene.")
@click.option("--uniform", is_flag=True, help="Place B genes uniformly instead of in islands.")
@click.option(
    "--min-frac",
    default=0.8,
    show_default=True,
    type=float,
    help="Call genes with at least this fraction on the read (0 = whole genes only).",
)
@click.option("--seed", default=20261114, show_default=True, type=int)
@click.option(
    "--out-prefix", required=True, help="Writes PREFIX.lengths.tsv.gz, PREFIX.features.tsv.gz, PREFIX.truth.tsv."
)
def simulate(
    read_lengths_path,
    a_lengths,
    b_lengths,
    a_class,
    b_class,
    n_genomes,
    total_bases,
    linked,
    uniform,
    min_frac,
    seed,
    out_prefix,
):
    """Simulate reads from a community with known linkage (truth computed on the genomes)."""
    import numpy as np

    from .io import read_lengths
    from .simulate import simulate_tables

    rl = read_lengths(read_lengths_path).length.to_numpy()
    rng = np.random.default_rng(seed + 999)

    def gene_lengths(path, median, sigma):
        if path:
            return np.loadtxt(path, dtype=np.int64, ndmin=1)
        return np.maximum(150, rng.lognormal(np.log(median), sigma, 5_000)).astype(np.int64)

    sim = simulate_tables(
        rl,
        gene_lengths(a_lengths, 900, 0.35),
        gene_lengths(b_lengths, 1000, 0.5),
        a_class=a_class,
        b_class=b_class,
        seed=seed,
        n_genomes=n_genomes,
        total_bases=total_bases,
        linked=linked,
        clustered=not uniform,
        min_frac=min_frac or None,
    )
    sim.lengths.to_csv(f"{out_prefix}.lengths.tsv.gz", sep="\t", index=False)
    sim.features.to_csv(f"{out_prefix}.features.tsv.gz", sep="\t", index=False)
    sim.truth.to_csv(f"{out_prefix}.truth.tsv", sep="\t", index=False)
    net = sim.truth[(sim.truth["mode"] == "end") & (sim.truth.w == 2000)].net.iloc[0]
    click.echo(
        f"{len(sim.lengths):,} reads, {len(sim.features):,} features; "
        f"true net linked (far-end distance, w = 2 kb): {net:.4f}"
    )


@main.command()
@click.option(
    "--lengths",
    "lengths_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Lengths table of a pilot dataset.",
)
@click.option(
    "--features",
    "features_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Feature table of the pilot dataset.",
)
@click.option("--a", "a_class", required=True, help="Class A (e.g. ARG).")
@click.option("--b", "b_class", required=True, help="Class B (e.g. MGE).")
@click.option("--bases", multiple=True, type=float, help="Depths to evaluate (default: 0.5, 1, 2, 5, 10, 20 Gbp).")
@click.option("--w", default=2000, show_default=True, type=int, help="Window of the net linked fraction (bp).")
@click.option("--net", default=0.1, show_default=True, type=float, help="Assumed net linked fraction.")
@click.option("--deff", default=1.0, show_default=True, type=float, help="Design effect for the half-width.")
@click.option("--target-halfwidth", default=0.05, show_default=True, type=float)
@click.option("--min-frac", default=0.8, show_default=True, type=float)
def power(lengths_path, features_path, a_class, b_class, bases, w, net, deff, target_halfwidth, min_frac):
    """What a given depth can detect, with densities and lengths taken from a pilot dataset (same read-length
    distribution assumed): detectable O/E of the read-level test and the half-width of the net linked fraction."""
    from .features import drop_same_orf
    from .io import read_features, read_lengths
    from .nulls.analytic import exposure, fit_density
    from .power import calls_density, depth_table, required_bases

    ln = read_lengths(lengths_path)
    f = read_features(features_path)
    fa, fb, _ = drop_same_orf(f[f["class"] == a_class], f[f["class"] == b_class])
    if fa.empty or fb.empty:
        raise click.ClickException("one of the classes has no features")
    L = ln.length.to_numpy()
    idx = ln.molecule_id
    has_a, has_b = idx.isin(fa.molecule_id).to_numpy(), idx.isin(fb.molecule_id).to_numpy()
    la, lb = (fa.end - fa.start).to_numpy(float), (fb.end - fb.start).to_numpy(float)
    lam_a = fit_density(exposure(L, la, min_frac), int(has_a.sum()))
    lam_b = fit_density(exposure(L, lb, min_frac), int(has_b.sum()))
    lam_a_calls, lam_b_calls = calls_density(len(fa), L, la), calls_density(len(fb), L, lb)
    pilot = float(L.sum())
    grid = sorted(set(bases or (5e8, 1e9, 2e9, 5e9, 1e10, 2e10)) | {pilot})
    t = depth_table(L, la, lb, lam_a, lam_b, lam_a_calls, lam_b_calls, grid, w=w, net=net, deff=deff, min_frac=min_frac)
    t.insert(0, "pilot", t.bases == pilot)
    click.echo(t.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    g = required_bases(L, la, lam_a_calls, lam_b_calls, w, target_halfwidth, net, deff)
    click.echo(f"bases for a net linked half-width of {target_halfwidth} at w = {w}: {g:.3g}")


@main.command("oe")
@click.option(
    "--lengths",
    "lengths_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Lengths table (molecule_id or read_id, length[, stratum column]).",
)
@click.option(
    "--features",
    "features_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Feature table (molecule_id, start, end, class[, scov]).",
)
@click.option("--pairs", required=True, multiple=True, help="Class pairs A:B, e.g. --pairs ARG:MGE --pairs ARG:MRG.")
@click.option("--strata-col", default=None, help="Column of the lengths table used as stratum (sample, taxon, ...).")
@click.option("--w", default=2000, show_default=True, type=int, help="Window of the net linked fraction (bp).")
@click.option("--boot", default=200, show_default=True, type=int, help="Bootstrap replicates (net linked fraction).")
@click.option(
    "--fragment",
    "targets",
    multiple=True,
    type=int,
    help="Also estimate the net linked fraction after random fragmentation to this N50 (repeatable).",
)
@click.option("--draws", default=5, show_default=True, type=int, help="Breakpoint draws per fragmentation target.")
@click.option(
    "--min-frac",
    default=0.8,
    show_default=True,
    type=float,
    help="Minimum fraction of a gene on the molecule for a call (exposure model).",
)
@click.option("--seed", default=20261114, show_default=True, type=int)
@click.option("--out", "out_path", required=True, help="Result table (TSV); settings go to OUT.meta.json.")
def oe_cmd(lengths_path, features_path, pairs, strata_col, w, boot, targets, draws, min_frac, seed, out_path):
    """Descriptive v0.1 table: molecule-level O/E (N2-ev, N0-ev, length-only) and the unstratified net linked fraction.

    Not the registered test (that is `linkstat test`); kept for comparison and for fragmentation (--fragment)."""
    import hashlib
    import json
    from datetime import datetime, timezone

    from .analysis import analyse_pairs
    from .io import read_features, read_lengths, validate_features

    ln = read_lengths(lengths_path)
    f = read_features(features_path)
    f = f[f.molecule_id.isin(ln.molecule_id)]
    validate_features(f, ln)
    if strata_col and strata_col not in ln.columns:
        raise click.UsageError(f"--strata-col {strata_col!r} is not a column of the lengths table")
    pair_list = []
    for p in pairs:
        a, sep, b = p.partition(":")
        if not sep or not a or not b:
            raise click.UsageError(f"pair {p!r} is not of the form A:B")
        pair_list.append((a, b))
    t = analyse_pairs(
        ln,
        f,
        pair_list,
        strata_col=strata_col,
        w=w,
        n_boot=boot,
        targets=tuple(targets),
        draws=draws,
        seed=seed,
        min_frac=min_frac,
    )
    t.to_csv(out_path, sep="\t", index=False)

    def sha(path):
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()

    meta = dict(
        linkstat=__version__,
        created=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        lengths=lengths_path,
        lengths_sha256=sha(lengths_path),
        features=features_path,
        features_sha256=sha(features_path),
        pairs=list(pairs),
        strata_col=strata_col,
        w=w,
        boot=boot,
        fragment=list(targets),
        draws=draws,
        min_frac=min_frac,
        seed=seed,
        null_model=(
            "N2-ev (exact permutation within bins of the excluded-volume exposure;"
            " N0-ev and the registered length-only null reported)"
        ),
        note="p values ignore reads that cover the same locus (pseudo-replication)",
    )
    with open(f"{out_path}.meta.json", "w") as fh:
        json.dump(meta, fh, indent=2)
    show = [
        "class_a",
        "class_b",
        "observed",
        "expected",
        "oe_ratio",
        "ci_low",
        "ci_high",
        "p",
        "q",
        "intact_net",
        "intact_lo",
        "intact_hi",
        "flags",
    ]
    click.echo(t[[c for c in show if c in t]].to_string(index=False, float_format=lambda x: f"{x:.4g}"))


@main.command("test")
@click.option("--lengths", "lengths_path", required=True, type=click.Path(exists=True, dir_okay=False),
              help="Lengths table (molecule_id or read_id, length[, strata column]).")
@click.option("--features", "features_path", required=True, type=click.Path(exists=True, dir_okay=False),
              help="Feature table (molecule_id, start, end, class[, scov]), e.g. from `linkstat prep`.")
@click.option("--pairs", required=True, multiple=True, help="Class pairs A:B, e.g. --pairs ARG:MGE --pairs ARG:MRG.")
@click.option("--kraken", "kraken_path", type=click.Path(exists=True, dir_okay=False),
              help="Kraken2 labels per read (read_id, taxid[, hit_kmers]), from all reads at confidence 0.")
@click.option("--taxonomy", "taxonomy_path", type=click.Path(exists=True, dir_okay=False),
              help="Kraken2 ktaxonomy.tsv of the database used (needed with --kraken).")
@click.option("--rank", default="S", show_default=True, type=click.Choice(["S", "G"]), help="Host stratum rank.")
@click.option("--min-hits", default=10, show_default=True, type=int,
              help="Reads with fewer Kraken2 hit k-mers are 'unresolved'.")
@click.option("--strata-col", default=None, help="Use this column of the lengths table as host stratum instead.")
@click.option("--loci", "paf_path", type=click.Path(exists=True, dir_okay=False),
              help="minimap2 -x ava-ont PAF among the carrier reads (unit-locus clumps).")
@click.option("--exclude", "exclude_path", type=click.Path(exists=True, dir_okay=False),
              help="Molecule ids to drop before the analysis: one per line, a table with a molecule_id/read_id "
                   "column, or a `linkstat chimeras` table (its flagged reads).")
@click.option("--w", default=2000, show_default=True, type=int, help="Window of net_w (bp).")
@click.option("--units-d", default=1000, show_default=True, type=int, help="A calls within d form one unit (bp).")
@click.option("--n-mc", default=9999, show_default=True, type=int, help="Monte Carlo draws of the null.")
@click.option("--no-oe", is_flag=True, help="Skip the descriptive unit-level O/E columns.")
@click.option("--seed", default=20261114, show_default=True, type=int,
              help="Seed for anchors and Monte Carlo draws (vary it per dataset in simulation studies).")
@click.option("--out", "out_path", required=True, help="Result table (TSV); settings go to OUT.meta.json.")
def test_cmd(lengths_path, features_path, pairs, kraken_path, taxonomy_path, rank, min_hits, strata_col, paf_path,
             exclude_path, w, units_d, n_mc, no_oe, seed, out_path):
    """The registered net linked test (amendment 04): is A linked to B beyond chance, and how strongly?

    One row per pair: net_w (within w) and net_mol (same molecule) with intervals, the two registered p values
    (Monte Carlo mid-p of the null-variance clump tests), BH q values across the pairs, flags, and the descriptive
    unit-level O/E."""
    import json
    from datetime import datetime, timezone

    from .hosts import read_taxonomy, strata_for_molecules
    from .io import read_features, read_lengths, validate_features
    from .loci import read_paf
    from .pairs import net_test_pairs

    ln = read_lengths(lengths_path)
    f = read_features(features_path)
    f = f[f.molecule_id.isin(ln.molecule_id)]
    validate_features(f, ln)
    pair_list = []
    for p in pairs:
        a, sep, b = p.partition(":")
        if not sep or not a or not b:
            raise click.UsageError(f"pair {p!r} is not of the form A:B")
        pair_list.append((a, b))
    strata, strata_src = None, None
    if kraken_path and strata_col:
        raise click.UsageError("use either --kraken or --strata-col")
    if kraken_path:
        if not taxonomy_path:
            raise click.UsageError("--kraken needs --taxonomy (ktaxonomy.tsv of the Kraken2 database)")
        k = pd.read_csv(kraken_path, sep="\t", dtype={"read_id": str})
        strata = strata_for_molecules(ln.molecule_id.to_numpy(), k, read_taxonomy(taxonomy_path), rank=rank,
                                      min_hits=min_hits)
        strata_src = f"kraken2 {kraken_path} rank {rank} min_hits {min_hits}"
    elif strata_col:
        if strata_col not in ln.columns:
            raise click.UsageError(f"--strata-col {strata_col!r} is not a column of the lengths table")
        strata, strata_src = ln[strata_col].to_numpy(), f"column {strata_col}"
    paf = read_paf(paf_path, min_block=1000) if paf_path else None
    exclude = None
    if exclude_path:
        e = pd.read_csv(exclude_path, sep="\t", header=None, dtype=str)
        col = 0
        if e.iloc[0].isin(["molecule_id", "read_id"]).any():
            head = e.iloc[0].tolist()
            e = e.iloc[1:]
            col = head.index("molecule_id") if "molecule_id" in head else head.index("read_id")
            if "flag" in head:  # a `linkstat chimeras` table: drop the flagged reads only
                e = e[e[head.index("flag")].str.lower().isin(["true", "1"])]
        exclude = set(e[col])
    t = net_test_pairs(ln, f, pair_list, strata=strata, paf=paf, exclude=exclude, w=w, units_d=units_d, seed=seed,
                       n_mc=n_mc, oe=not no_oe)
    t.to_csv(out_path, sep="\t", index=False)
    meta = dict(
        linkstat=__version__,
        created=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        lengths=lengths_path, lengths_sha256=_sha256(lengths_path),
        features=features_path, features_sha256=_sha256(features_path),
        pairs=list(pairs), strata=strata_src, loci=paf_path, exclude=exclude_path,
        excluded_molecules=len(exclude) if exclude else 0,
        w=w, units_d=units_d, n_mc=n_mc, seed=seed,
        test="net linked test (amendment 04): null-variance clump score tests, Monte Carlo mid-p; BH per test",
    )
    with open(f"{out_path}.meta.json", "w") as fh:
        json.dump(meta, fh, indent=2)
    show = ["class_a", "class_b", "net_w", "net_w_lo", "net_w_hi", "net_w_p", "net_w_q", "net_mol", "net_mol_lo",
            "net_mol_hi", "net_mol_p", "net_mol_q", "flags"]
    click.echo(t[[c for c in show if c in t]].to_string(index=False, float_format=lambda x: f"{x:.4g}"))


@main.command("chimeras")
@click.option("--paf", "paf_path", type=click.Path(exists=True, dir_okay=False),
              help="minimap2 -x ava-ont PAF among the carrier reads (overlap breakpoints).")
@click.option("--reads", "reads_path", type=click.Path(exists=True, dir_okay=False),
              help="Reads (FASTA/FASTQ, gzip ok) to scan for adapters inside reads (needs edlib).")
@click.option("--barcodes", "barcodes_path", type=click.Path(exists=True, dir_okay=False),
              help="Native barcodes to scan for, as a TSV of name and sequence (e.g. NB01 CACAAAG...).")
@click.option("--min-support", default=2, show_default=True, type=int, help="Partner reads supporting a breakpoint.")
@click.option("--min-edge", default=500, show_default=True, type=int, help="Breakpoint distance from read ends (bp).")
@click.option("--out", "out_path", required=True, help="Flag table (TSV) for `linkstat test --exclude`.")
def chimeras_cmd(paf_path, reads_path, barcodes_path, min_support, min_edge, out_path):
    """Flag likely chimeric reads without a reference: overlap breakpoints and adapters inside reads.

    Mock communities (P2 step 2.5): the defaults flag 0.6-2.2% of clean reads and 73-100% of reads joining two species.
    The registered test does not drop them by itself; use the table with `linkstat test --exclude`."""
    from .chimeras import chimera_flags, read_sequences
    from .loci import read_paf

    if not paf_path and not reads_path:
        raise click.UsageError("give --paf, --reads, or both")
    bc = None
    if barcodes_path:
        b = pd.read_csv(barcodes_path, sep="\t", header=None, names=["name", "seq"], dtype=str)
        bc = dict(zip(b.name, b.seq, strict=True))
    paf = read_paf(paf_path, min_block=1000) if paf_path else None
    reads = read_sequences(reads_path) if reads_path else None
    t = chimera_flags(paf=paf, reads=reads, barcodes=bc, min_support=min_support, min_edge=min_edge)
    t.to_csv(out_path, sep="\t", index=False)
    click.echo(f"{len(t)} reads judged, {int(t.flag.sum())} flagged "
               f"({int(t.flag_overlap.sum())} overlap breakpoints, {int(t.flag_adapter.sum())} adapters inside)")


def _sha256(path):
    import hashlib

    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


if __name__ == "__main__":
    main()
