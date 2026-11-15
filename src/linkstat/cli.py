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


if __name__ == "__main__":
    main()
