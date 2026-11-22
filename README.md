# mobilink

mobilink tests whether two gene classes sit close together on long reads or contigs more often than chance, and
measures how strongly they are linked with a value that does not depend on read length.

It was built for antibiotic-resistance genes (ARG), metal-resistance genes (MRG) and mobile genetic elements (MGE) in
nanopore metagenomes, but it accepts any two classes of features on any molecules.

**Status: release candidate (v1.0.0.dev0).** Not yet on PyPI or Bioconda.

## Why

The usual measure is the share of ARG-carrying reads that also carry an MGE. That share grows with read length even
when the genes are placed at random, so it cannot be compared between studies. mobilink replaces it with two things:

- **The net linked fraction**, an effect size. It is the share of A genes that have a B gene within a distance *w*
  (default 2 kb), beyond what chance placement gives. It is 0 under chance and stays the same when reads get shorter.
- **A test against chance** that keeps the A genes where they are, respects the host of each read, and counts reads
  covering the same genomic locus as one observation.

In 14 public runs (9 human gut, 4 livestock, 1 wastewater), cutting the reads to 3 kb lowered the naive share 1.8- to
7.7-fold, while the net linked fraction changed by at most 0.006 (see `benchmark/`).

## Install

```bash
pip install .            # from a clone of this repository; Python 3.10 or later
mobilink --version
```

The only dependencies are numpy, scipy, pandas and click. `mobilink chimeras --reads` also needs edlib
(`pip install ".[adapters]"`).

## Five-minute example

The repository holds a read-length table (`examples/read_lengths.tsv.gz`, 20,000 lengths from a public human gut
run, read N50 18.5 kb). The commands below simulate a community in which 30% of the ARG sit next to an MGE, run the
test, and then cut the reads.

```bash
mobilink simulate --read-lengths examples/read_lengths.tsv.gz --n-genomes 60 --total-bases 1e9 --linked 0.3 \
  --out-prefix toy
mobilink test --lengths toy.lengths.tsv.gz --features toy.features.tsv.gz --pairs ARG:MGE --strata-col genome \
  --out toy.results.tsv
mobilink oe --lengths toy.lengths.tsv.gz --features toy.features.tsv.gz --pairs ARG:MGE --strata-col genome \
  --fragment 9000 --fragment 4000 --out toy.cut.tsv
```

| Read N50 | Naive share | Net linked fraction within 2 kb |
|---|---|---|
| 18.5 kb (as simulated) | 0.42 | 0.131 |
| cut to 9 kb | 0.31 | 0.130 |
| cut to 4 kb | 0.20 | 0.135 |

The true value of the simulated community is 0.125 (`toy.truth.tsv`, row `end`, w = 2000). `docs/worked_example.md`
explains every step and every number.

## Real data

mobilink needs two tab-separated tables:

- **lengths:** `molecule_id` (or `read_id`), `length`, and optionally a host column;
- **features:** `molecule_id`, `start`, `end`, `class` (0-based, end exclusive), and optionally `scov` (subject cover
  in %).

`mobilink prep` builds the feature table from DIAMOND blastx hits against CARD, BacMet2 and mobileOG-db. For reads,
the test also takes Kraken2 labels (host strata) and minimap2 overlaps among the gene-carrying reads (locus clumps):

```bash
mobilink test --lengths reads.lengths.tsv.gz --features reads.primary.features.tsv.gz \
  --pairs ARG:MGE --pairs ARG:MRG --pairs MRG:MGE \
  --kraken reads.kraken.tsv.gz --taxonomy k2db/ktaxonomy.tsv --loci carriers.ava.paf.gz --out results.tsv
```

For contigs, give the bin of each contig as the host column and use `--molecule-clumps`:

```bash
mobilink test --lengths contigs.lengths.tsv --features contigs.features.tsv --pairs ARG:MGE \
  --strata-col bin --molecule-clumps --out contig_results.tsv
```

`docs/real_data.md` gives every command from the FASTQ file to the result table.

## Commands

| Command | Purpose |
|---|---|
| `mobilink prep` | Feature table from DIAMOND hits (one call per ORF and database) |
| `mobilink test` | The net linked test: effect sizes, intervals, p and q values, flags |
| `mobilink meta` | Estimates for groups of runs (e.g. environments) with a delete-one-project jackknife, and contrasts |
| `mobilink chimeras` | Flags likely chimeric reads without a reference, for `test --exclude` |
| `mobilink power` | What a given sequencing depth can detect, from a pilot dataset |
| `mobilink simulate` | Reads from a community with known linkage |
| `mobilink oe` | Older descriptive table; also cuts reads (`--fragment`) |

## Reading the result

The main columns of `mobilink test` are:

| Column | Meaning |
|---|---|
| `net_w`, `net_w_lo`, `net_w_hi` | Net linked fraction within *w*, with a 95% interval (delete-one-host jackknife) |
| `net_w_p`, `net_w_q` | Test of `net_w` against chance; q is Benjamini–Hochberg over the pairs |
| `net_w_O`, `net_w_E` | Observed and expected events within *w* |
| `net_mol`, `net_mol_p` | The same for "B anywhere on the same molecule" |
| `mol_ratio`, `mol_ratio_lo90`, `mol_ratio_hi90` | Share of A units with a B on the molecule, against chance, as a ratio |
| `hosts_effective` | Effective number of hosts behind the estimate |
| `flags` | Reasons to distrust the row |

Read the flags first. An empty estimate with the flag "no events" or "few effective events" means the data cannot
answer the question. It does not mean that the classes avoid each other. `docs/output_columns.md` lists all columns
and flags.

## Limits

- **Power.** Soil reads with an N50 near 5 kb give almost no ARG–MGE events: 10 Gbp gave 0.36 expected chance events
  in our 27 soil libraries. Run `mobilink power` on a pilot before sequencing.
- **Few hosts.** When one host dominates a sample, the intervals are wide. This is correct, not a fault.
- **Class definitions.** The value depends on what counts as A and B (identity thresholds, which mobile-element
  categories). Report the class definitions with every estimate.
- **Cut reads.** `mobilink oe --fragment` keeps a gene that is cut by a breakpoint when 80% of it remains. Long genes
  are then slightly over-represented, and in 8 public runs this lowered the estimate by 0.5–6% in relative terms.
  This belongs to the cut experiment, not to the estimate on native reads; with whole-gene retention the change was
  centred on zero.
- **Contigs.** "Same contig" depends on how fragmented the assembly is.

## Method, validation and citation

- `docs/method.md` describes the estimator and the test.
- `benchmark/` holds the scripts and tables behind the validation: simulations with known truth (5.0% and 4.4% false
  positives at a nominal 5% in 9,000 simulated datasets without linkage), mock communities, and the pre-registered tests on public
  runs.
- A paper is in preparation. Until it is out, cite this repository and the version you used.

## Development

```bash
pip install -e ".[dev]"
pytest
```

## Licence

MIT; see `LICENSE`. The helper tables in `resources/` are derived from mobileOG-db and BacMet2; see
`resources/README.md` for their sources.
