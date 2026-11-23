# Phase 4, step 4.9: package version 1.0 (candidate)

Date: 2026-11-22. State: developed on branch `v1.0-dev` and, after the data passes of steps 4.7 and 4.7b had ended,
merged into the main branch of the local package repository (`pkg/`, commit 0c86ffd, version 1.0.0.dev0; 114 of 114
tests pass on main). Not pushed, not tagged, not published. The analyses of the paper used commit 86e3696 (version
0.5.0.dev0).

## What changed

| Change | Why | Effect on results |
|---|---|---|
| `mobilink test`: the saddlepoint system of the descriptive O/E columns builds its matrix over the used coordinates only | With about 10,000 host strata it built a dense matrix over all strata: 9.7 GB and most of the run time | None. Result table byte-identical on the run checked; all 64 release-check runs exact |
| `mobilink chimeras`: the overlap file is read in one pass; small records go to temporary files (`--tmp-dir`) | All overlaps were held in memory: up to 104 GB | None from the reader itself |
| `mobilink chimeras`: fix of ties | When two overlaps of a read started at the same position, the veto of a breakpoint by a spanning overlap depended on the row order and could be missed | 1,984 of 674,080 flags removed in 64 runs (0.3%; at most 0.9% of a run; 62 runs affected); none added; the same reads are judged |
| New command `mobilink meta` | The model of the survey (pooled class estimate, delete-one-project jackknife, contrasts) existed only in the frozen analysis script | None (new) |
| Clear message when no class A call is left after the same-ORF rule | The case stopped the sub-class run with a numpy error text | None |

## Checks

- 114 of 114 tests pass; the linter is clean.
- **Release check** (`scripts/P4_09_release.sh`, `tables/P4/v1check/P4_09_release_check.tsv`): `mobilink test` of the
  branch against the registered per-run tables of the 64 Phase 3 runs: 192 pairs, 58 numeric columns, largest
  difference 0, missing values and flags equal in 64 of 64. No frozen analysis was rerun; the registered tables were
  only read.
- **Streamed overlaps** equal the in-memory result on a real run (ERR7625327, 114,758 judged reads) and in the
  tests (three chunk and bucket settings).
- **`meta` against the registered class estimates** (all 241 runs; `tables/P4/v1check/meta_ARG_MGE.*.tsv`): the seven
  class estimates with standard error, interval, degrees of freedom and effective number of projects, the 102
  jackknife replicates and the 21 contrasts (difference, standard error, degrees of freedom, t, p, Holm p) equal the
  registered tables exactly (largest difference 0). One process, 127 CPU minutes, 20.3 GB at the peak (the human gut
  class, 64 runs). The registered tables were only read.

## Run time and memory, 64 public runs, one thread

| Command | Version 0.5: CPU min | Version 0.5: peak GB | Version 1.0: CPU min | Version 1.0: peak GB |
|---|---|---|---|---|
| `prep` | 0.37 (0.09–11.3) | 1.2 (0.5–25.8) | 0.21 (0.05–6.3) | 1.2 (0.5–25.8) |
| `chimeras` | 0.21 (0.01–20.9) | 1.0 (0.1–104) | 0.08 (0.02–7.2) | 0.55 (0.08–1.1) |
| `test` | 3.0 (0.6–12.9) | 7.8 (0.5–14.2) | 0.81 (0.19–7.7) | 1.2 (0.34–2.1) |

Median (range). The two measurements were made under different server loads (about 100–130 and about 20–30), so the
CPU times are not a clean comparison: `prep` did not change and took about half the CPU time in the second
measurement. The memory figures do not depend on the load. Tables: `tables/P3/P3_06_runtime.tsv` and
`tables/P4/v1check/P4_09_runtime.tsv`.

## Open

- The release itself (tag, repository, PyPI, Bioconda, Zenodo) is the user's decision.
- `/tmp` on this server is held in memory; `mobilink chimeras` should be given `--tmp-dir` on disk here.
