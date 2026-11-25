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

This downloads CARD (the current release), BacMet2 (experimentally confirmed), mobileOG-db beatrix-1.6 (355 MB) and
NCBI UniVec (vector sequences, 2 MB), checks them, and builds the DIAMOND indexes: about 1 GB of disk and a few minutes. The versions, checksums, licences
and citations are written to `mobilink_db/mobilink_db.json` and repeated in every report.

- **Licences.** CARD is free for academic research; commercial use needs a licence from McMaster University.
  mobileOG-db is CC BY 4.0. Cite the three databases when you publish.
- **Without internet access.** Download the files elsewhere and give them with `--card` (the CARD data archive or
  the protein homolog FASTA), `--bacmet`, `--mobileog` and `--univec`.
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

Contigs are searched in windows of 30 kb, because a search program reports one alignment per reference protein and
query by default and would miss the second copy of a repeated gene, such as an insertion sequence, on a long contig.
The hit tables of a contig run are in contig coordinates.

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
| `NAME.power.tsv` | Per pair: the half-width of the interval, the smallest net linked fraction the test finds with 80% power at this depth, and the sequence needed to detect 0.05, 0.02 and 0.01 (reads only) |
| `NAME.fixed_length.tsv` | Screen for genes whose carrier reads share one read length, with the UniVec check of each flagged gene (reads only) |
| `NAME.test_no_lab_plasmid.tsv` | The result without the reads of a laboratory plasmid; written only when one is found |
| `NAME.features.tsv.gz`, `NAME.lengths.tsv.gz`, `NAME.kraken.tsv.gz`, `NAME.ava.paf.gz`, `NAME.chimeras.tsv` | The tables behind the result; they are the inputs of `mobilink test`, `mobilink meta` and `mobilink chimeras` |
| `NAME.run.json` | Settings, versions, database checksums, input numbers and flags |

## Reading the report

- **The net linked fraction** counts the two flanks of each gene unit. It is the share of flanks with a whole gene of
  the other class within the window (2 kb by default), beyond what chance gives. A unit with such a gene on one flank
  only counts one half. Under chance the value is 0 at any read length.
- **"No linkage shown" is not "no linkage".** The section "What this depth can show" gives the smallest net linked
  fraction that the test would have found at this depth, and the sequence needed to find weaker linkage. `mobilink
  power` does the same from a pilot before sequencing.
- **The interval and the test answer different questions.** The test asks whether the genes sit closer than random
  positions in the same hosts; its detection limit falls with the square root of the depth. The interval comes from
  leaving out one host at a time, so it also carries the differences between hosts. On eight public runs thinned
  to one eighth of their reads the interval widened only 1.2-fold at the median, while the detection limit followed
  the depth. More reads of one sample therefore make weaker linkage detectable, but narrow the interval little. A
  result with a small p value and an interval that includes 0 means: linked in some hosts, not shown for all.
- **Host labels.** Without `--kraken-db`, genes that are common in one abundant host can look linked. With it, the
  report gives the share of molecules without a label. Outside the human gut this share is often above 70%, and the
  host control is then weak.
- **Laboratory plasmids.** If at least 30 reads carry one ARG and a quarter of them have the same length (within 3%)
  below 20 kb, the gene is flagged: this is the pattern of a small plasmid sequenced whole. The four carrier reads
  nearest that length are then aligned to NCBI UniVec. If at the median at least 35% of the read matches and at least
  two different vector records match with 300 bp or more, the cluster is called a laboratory plasmid, and the report
  adds a result table without any read that carries the gene. Use that table. A flagged gene that fails the check is
  probably a natural small plasmid and stays in. The rule was fixed on the 14 flagged clusters of the survey, where
  one soil project held an expression plasmid that gave 44% of all linked events of its class.
- **Chimeras.** The main result keeps all reads, as in the paper. If the table without flagged reads differs
  clearly, report both.

## Several samples

Run each sample, then list the per-sample tables in a manifest and use `mobilink meta` to estimate groups (for
example environments) with a jackknife over projects. See `real_data.md`, section "Groups of runs".
