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


if __name__ == "__main__":
    main()
