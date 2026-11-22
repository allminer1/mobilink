# Changelog

## 1.0.0 (unreleased)

Changes since commit 86e3696 (0.5.0.dev0), the version used for the survey of 241 public runs (amendment 07).

- New command `mobilink meta` (module `mobilink.meta`): one estimate per group of runs, with a delete-one-cluster
  jackknife interval, the effective number of clusters, and contrasts between groups (Welch on the jackknife standard
  errors, Holm). It is the model of the survey in the paper. `mobilink.stats.nettest.net_w_estimate` gives the
  `net_w` of the test without the tests, for replicates.
- `mobilink chimeras` reads the overlap file in one pass and no longer holds it in memory (new option `--tmp-dir`;
  44 bytes of temporary disk space per usable overlap). A run with 445 million overlap rows needed 0.65 GB and
  7.5 CPU minutes in place of 104 GB and 21 CPU minutes.
- `mobilink chimeras`, fix: when two overlaps of a read started at the same position, the veto of a breakpoint by a
  spanning overlap depended on the order of the rows and could be missed. The result no longer depends on the row
  order. Effect on two public runs: 1 of 13,012 and 5 of 42,392 flags removed; none added.
- `mobilink test`: the descriptive O/E columns need far less memory and time when there are thousands of host strata
  (the saddlepoint system built a dense matrix over all strata). The result table is unchanged. A soil run with 0.98
  million reads and 9,883 strata: 0.7 GB and 82 s of CPU time in place of 12.9 GB and 152 s.
- A pair in which every class A call shares its open reading frame with a class B call now gives the flag "no class A
  call is left after the same-ORF rule" (before: an error text from numpy).

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
