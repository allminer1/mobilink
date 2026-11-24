# Helper tables for `mobilink prep` and `mobilink run`

The two tables are installed with the package (`src/mobilink/data/`), and `mobilink db` copies them into the database
directory. The two scripts that built them are here.

| File | What it is | Source |
|---|---|---|
| `mobility_core.tsv.gz` | The 109,325 mobileOG-db entries that count as MGE in the primary rule (the "mobility core") | mobileOG-db beatrix-1.6, built with `make_mobility_core.py` |
| `bacmet2_exp_classes.tsv` | The 753 BacMet2 experimentally confirmed entries, with the metals each one names and `is_mrg` (TRUE for 489 entries that name at least one metal or metalloid; the other 264 are biocide-only) | BacMet 2.0 mapping file `BacMet2_EXP.753.mapping.txt`, column "Compound" |
| `make_mobility_core.py` | Builds the mobility core from the mobileOG-db metadata table | |
| `make_mobileog_lookup.py` | Builds the lookup table for `--mobileog-lookup` from the mobileOG-db FASTA file | |

The sequence databases themselves are not included. `mobilink db` downloads them; cite them:

- CARD (protein homolog model): https://card.mcmaster.ca
- mobileOG-db (beatrix-1.6): https://mobileogdb.flsi.cloud.vt.edu (archived copy: Zenodo record 13241605)
- BacMet 2.0 (experimentally confirmed): http://bacmet.biomedicine.gu.se

## The mobility core

mobileOG-db groups proteins of mobile elements into broad categories, and many of them have chromosomal homologs
(recombination and repair, replication, partitioning, regulators). Counting all of them as MGE marks a large share of
ordinary chromosomal reads as "mobile". The mobility core keeps the entries that move DNA within or between cells:

1. major category integration/excision or transfer;
2. evidence Manual or Homology (not Keyword Search);
3. no minor category among regulation, partitioning, competence, replication/recombination/repair;
4. gene family not in a fixed list of 42 chromosomal housekeeping families (listed in `make_mobility_core.py`).

The rule uses only the annotation fields of the database, never the hits of a dataset.
