# From a FASTQ file to a report in two commands

`mobilink run` does the whole analysis of one sample: gene calls, host labels, overlaps, the net linked test, a
per-gene table and an HTML report. It uses the settings of the survey in the paper. The step-by-step commands that
it replaces are in `real_data.md`; the estimator and the test are those of `mobilink test`.

## What you need

| Program | Used for | Needed |
|---|---|---|
| DIAMOND 2.1 or later | Gene calls (blastx, long-read mode) | always |
| minimap2 | Overlaps among the reads that carry a gene | for reads |
| Kraken2 and a prebuilt index | Host labels | recommended |

```bash
mamba install -c conda-forge -c bioconda mobilink diamond minimap2 kraken2
```

## 1. Databases (once)

```bash
mobilink db --out mobilink_db --threads 8
```

This downloads CARD (the current release), BacMet2 (experimentally confirmed) and mobileOG-db beatrix-1.6 (355 MB),
checks them, and builds the DIAMOND indexes: about 1 GB of disk and a few minutes. The versions, checksums, licences
and citations are written to `mobilink_db/mobilink_db.json` and repeated in every report.

- **Licences.** CARD is free for academic research; commercial use needs a licence from McMaster University.
  mobileOG-db is CC BY 4.0. Cite the three databases when you publish.
- **Without internet access.** Download the files elsewhere and give them with `--card` (the CARD data archive or
  the protein homolog FASTA), `--bacmet` and `--mobileog`.
- **Host labels.** `mobilink db` does not fetch a Kraken2 index. Download a prebuilt one from the Kraken2 index page
  (the paper used Standard-16, 16 GB) and unpack it; the directory must hold `ktaxonomy.tsv`.

## 2. One sample

```bash
mobilink run reads.fastq.gz --db mobilink_db --kraken-db k2_standard_16gb --out results --threads 16
```

For an assembly, add `--contigs` (no overlap and chimera steps; all gene units of one contig count as one locus):

```bash
mobilink run contigs.fa --contigs --db mobilink_db --kraken-db k2_standard_16gb --out results_contigs
```

A stopped run resumes where it stopped: a step whose output file exists is skipped. To redo a step, delete its file.

## What you get

`NAME` is the input file name without its extension, or `--name`.

| File | Content |
|---|---|
| `NAME.report.html` | The report: one sentence per pair, the result table, the distance curve, input numbers, flags, the genes behind the result. One file, no internet needed to open it |
| `NAME.test.tsv` | The result table of `mobilink test` (columns: `output_columns.md`) |
| `NAME.test_no_chimeras.tsv` | The same without reads flagged as chimeras (reads only) |
| `NAME.curve.tsv` | Net linked fraction for windows of 0.25 to 10 kb |
| `NAME.genes.tsv` | One row per gene: calls, calls with a gene of the other class within the window, the commonest neighbour, and for ARG the CARD mechanism and drug class |
| `NAME.fixed_length.tsv` | Screen for genes whose carrier reads share one read length (reads only) |
| `NAME.features.tsv.gz`, `NAME.lengths.tsv.gz`, `NAME.kraken.tsv.gz`, `NAME.ava.paf.gz`, `NAME.chimeras.tsv` | The tables behind the result; they are the inputs of `mobilink test`, `mobilink meta` and `mobilink chimeras` |
| `NAME.run.json` | Settings, versions, database checksums, input numbers and flags |

## Reading the report

- **The net linked fraction** counts the two flanks of each gene unit. It is the share of flanks with a whole gene of
  the other class within the window (2 kb by default), beyond what chance gives. A unit with such a gene on one flank
  only counts one half. Under chance the value is 0 at any read length.
- **"No linkage shown" is not "no linkage".** Look at the interval. A wide interval means the sample is too shallow
  for this pair; `mobilink power` says how much sequence is needed.
- **Host labels.** Without `--kraken-db`, genes that are common in one abundant host can look linked. With it, the
  report gives the share of molecules without a label. Outside the human gut this share is often above 70%, and the
  host control is then weak.
- **Fixed-length flag.** If at least 30 reads carry one ARG and a quarter of them have the same length (within 3%)
  below 20 kb, the report says so. This is the pattern of a small plasmid sequenced whole. In one soil project of
  the survey it was a laboratory expression plasmid that gave 44% of all linked events of its class. Align a few of
  the flagged reads to NCBI UniVec before you interpret the ARG results.
- **Chimeras.** The main result keeps all reads, as in the paper. If the table without flagged reads differs
  clearly, report both.

## Several samples

Run each sample, then list the per-sample tables in a manifest and use `mobilink meta` to estimate groups (for
example environments) with a jackknife over projects. See `real_data.md`, section "Groups of runs".
