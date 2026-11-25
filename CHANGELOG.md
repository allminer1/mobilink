# Changelog

## 1.1.0 (in development, branch `wp3-run`)

The estimator and the test are unchanged from 1.0.0.

- New command `mobilink db`: downloads CARD (current release), BacMet2 EXP and mobileOG-db beatrix-1.6 (archived copy
  on Zenodo, checksum verified), builds the DIAMOND indexes and writes a manifest with versions, checksums, licences
  and citations. Local files can be given instead of downloads.
- New command `mobilink run`: one sample from reads or contigs to the net linked test, with host labels (Kraken2),
  overlaps (minimap2), chimera flags, a distance curve, a per-gene table (`mobilink.genes.gene_table`), the screen
  for fixed-length carrier reads (`mobilink.genes.fixed_length_clusters`) with a check of each flagged gene against
  NCBI UniVec (`mobilink.genes.vector_rule`; a laboratory plasmid gives a second result table without its reads),
  and a self-contained HTML report. A stopped run resumes.
- The report and `NAME.power.tsv` say what the depth of the run can show (`mobilink.power.depth_projection`): the
  half-width of the interval, the smallest net linked fraction the test finds with 80% power, and the sequence
  needed to detect weaker linkage. The depth columns use the spread of the test and not the interval: on eight
  public runs thinned to 1/8 the test variance followed the depth, the delete-one-host interval did not.
- `mobilink run --contigs` searches contigs in windows of 30 kb (every 20 kb) and keeps every alignment of a
  reference protein. Searched whole, long contigs lost the second and later copies of repeated genes: on 21 complete
  genomes, 685 mobility-core calls in place of 967. Found with the independent simulation, before any release.
- `mobilink run` keeps every alignment of a reference protein for reads as well (`--max-hsps 0`), so that a second
  copy of a repeated gene on one read is called. `--one-alignment` restores the calling of the analyses of the
  paper. Measured on the same reads of two public human gut runs (0.4 Gbp each, nanopore and PacBio HiFi): 1.0% and
  0.4% more calls of the mobility core, 0.3% and 0% more ARG calls, and an ARG–MGE net linked fraction that moved
  from 0.0483 to 0.0485 and from 0.0649 to 0.0667.
- `mobilink db` resumes a download that the server or the line cuts, and checks the size against the size announced.
  Before, a cut download of mobileOG-db (355 MB) was only found by its checksum, and the command had to be started
  again. Found by running the tutorial from the start on a slow line.
- The depth table (`NAME.power.tsv`, report section "What this depth can show") gives the half-width of the interval
  that the result table reports. Before any release it used 1.96 standard errors, which is narrower than the reported
  interval when few hosts carry the genes.
- Fixed before any release: a built wheel left out the mobility core, because the ignore rule for `*.tsv.gz` also
  applied to `src/mobilink/data/`. Editable installs were not affected, so the tests did not see it; the container
  check did. The continuous integration now checks the wheel.
- A `Dockerfile` (mobilink with DIAMOND, minimap2 and Kraken2; built and checked in the continuous integration) and
  a documentation site configuration (`mkdocs.yml`, `docs/index.md`).
- The BacMet2 class table and the mobility core are installed with the package (moved from `resources/` to
  `src/mobilink/data/`).
- `docs/output_columns.md`: `net_w` is described as a share of the sides of A units, as in `docs/method.md`.

## 1.0.0

Changes since commit 86e3696 (0.5.0.dev0), the version used for the survey of 241 public runs (amendment 07).
Commit numbers in this file are those cited in the registered files; each has a tag `orig-<number>`, and
`benchmark/COMMITS.tsv` gives the present numbers (see `benchmark/README.md`, "Commit numbers").

- New command `mobilink meta` (module `mobilink.meta`): one estimate per group of runs, with a delete-one-cluster
  jackknife interval, the effective number of clusters, and contrasts between groups (Welch on the jackknife standard
  errors, Holm). It is the model of the survey in the paper. `mobilink.stats.nettest.net_w_estimate` gives the
  `net_w` of the test without the tests, for replicates.
- `mobilink chimeras` reads the overlap file in one pass and no longer holds it in memory (new option `--tmp-dir`;
  44 bytes of temporary disk space per usable overlap). A run with 445 million overlap rows needed 0.65 GB and
  7.5 CPU minutes in place of 104 GB and 21 CPU minutes.
- `mobilink chimeras`, fix: when two overlaps of a read started at the same position, the veto of a breakpoint by a
  spanning overlap depended on the order of the rows and could be missed. The result no longer depends on the row
  order. Effect on 64 public runs: 1,984 of 674,080 flags removed (0.3%; at most 0.9% of the flags of a run), none
  added; the same reads are judged.
- `mobilink test`: the descriptive O/E columns need far less memory and time when there are thousands of host strata
  (the saddlepoint system built a dense matrix over all strata). The result table is unchanged. A soil run with 0.98
  million reads and 9,883 strata: 0.7 GB and 82 s of CPU time in place of 12.9 GB and 152 s.
- A pair in which every class A call shares its open reading frame with a class B call now gives the flag "no class A
  call is left after the same-ORF rule" (before: an error text from numpy).
- Release check: `mobilink test` reproduces the registered per-run tables of 64 public runs exactly (192 pairs, 58
  numeric columns, largest difference 0, flags equal).

## 0.5.0 (commit 86e3696)

- The package is renamed from the working name `linkstat` to `mobilink` (command, import name, and the version key
  in the `.meta.json` files).
- `mobilink test --molecule-clumps` (API: `loci="molecule"`): all A units of one molecule are one clump. Use it for
  contigs.
- `mobilink test --loci`: only overlaps among molecules with a class-A call are loaded. Results do not change; a
  large overlap file needs much less memory.
- `mobilink power --bases`: the help states the unit (bases).
- Documentation: worked example, real-data guide, output columns, method note; helper tables in `resources/`.
- Benchmark layout (`benchmark/`), CI workflow, Bioconda recipe draft.

## 0.3.0.dev0 (commit 9819182)

The version used for the pre-registered analyses (amendments 05 and 06). Same-molecule ratio with a 90% interval,
pooled runs, filtered overlap reader.

## 0.2.0.dev0 (commit cf0db83)

`test` runs the registered net linked test (amendment 04); chimera flags; reporting fixes.

## 0.1

Molecule-level O/E tables and the first net linked estimator.
