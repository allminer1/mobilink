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


@main.command("test")
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
def test_cmd(lengths_path, features_path, pairs, strata_col, w, boot, targets, draws, min_frac, seed, out_path):
    """Test whether class pairs share molecules more often than chance, and estimate how strongly they are linked."""
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
        null_model="N0-ev (excluded-volume corrected exposure null; registered length-only null reported)",
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


if __name__ == "__main__":
    main()
