# P3 steps 3.5 and 3.6: release package, run time and memory

2026-11-18. The package is named **mobilink** (user decision D85). Everything here is local: nothing is pushed or
published. Decisions D85–D94 in `docs/analysis_plan.md`.

## Result

**The release code gives exactly the registered results, and it runs in minutes per sample.**

- **Release check.** `mobilink test` (commit `a6b44c4` and later) reproduces the registered per-run table of every
  accepted public run of both stages: 64 of 64 runs, 3 pairs each, 58 numeric columns, largest difference 0, flags
  identical (`tables/P3/P3_06_release_check.tsv`). `mobilink prep` gives the same ARG, MRG and MGE calls as the frozen
  pipeline (checked on SRR29980985: 1,564 / 15 / 25,779 calls, identical coordinates).
- **Tests.** 104 of 104 pass (`logs/P3/pkg_tests_v05.log`). The lint step of the CI file could not be run here
  (`ruff` is not installed).
- **No frozen analysis was rerun.** The measurements write to `tmp/p35/`; nothing under `tables/P3/confirm*` or
  `public/` changed.

## Run time and memory (step 3.5, V8)

64 public runs of about 3 Gbp (0.24–1.9 million reads; 1.1–93 million DIAMOND hit rows; 0.02–445 million overlap
rows). One thread per command, `nice 15`, on a server with a load of 100–130 on 64 cores, so only CPU time (user +
system) and peak memory are reported. Tables: `tables/P3/P3_06_runtime.tsv` (per run and command) and
`P3_06_runtime.summary.tsv`.

| Command | CPU minutes, median (range) | Peak memory in GB, median (range) | CPU seconds per million reads, median | What it grows with |
|---|---|---|---|---|
| `prep` | 0.37 (0.09–11.3) | 1.2 (0.5–25.8) | 33 | DIAMOND hit rows: 0.30 GB and 5.2 s per million rows |
| `chimeras` | 0.21 (0.01–20.9) | 1.0 (0.1–104) | 19 | overlap rows: 0.34 GB and 4.1 s per million rows |
| `test` | 3.0 (0.6–12.9) | 7.8 (0.5–14.2) | 271 | host strata, through the descriptive O/E columns |

- **`test` is dominated by the descriptive O/E columns.** With `--no-oe`, soil run ERR11593393 (984,000 reads) needs
  0.5 GB and 11 s of user time in place of 12.9 GB and 223 s. The estimates, intervals and p values are identical
  (49 shared numeric columns, difference 0; also checked on SRR23278716). The default keeps the columns, because the
  registered per-run tables contain them.
- **`chimeras` loads every overlap.** 53 of 64 runs need less than 4 GB. Two runs need 85 GB (SRR33225186, 188
  million overlap rows) and 104 GB (SRR23278728, 445 million rows). This is a limit of the release; a streaming
  reader would remove it. The guide states it.
- **`prep`** needs 25.8 GB for one run with 93 million hit rows (SRR23278728); all others stay under 3 GB.
- The annotation before mobilink costs far more: DIAMOND alone ran at about 32 Gbp per hour on 16 threads per
  database (step 0.2).

## Three exploratory checks (outside the registration)

| Question | Answer | Table |
|---|---|---|
| Why is the cut estimate 0.5–6% lower in the 8 counting runs with a net linked fraction ≥ 0.06? | The cut keeps a gene when 80% of it remains, but the weights assume whole genes. Long ARG are then over-represented, and they are less often linked to MGE (7 of 8 runs). With whole-gene retention the relative change is centred on zero (mean 0.0001; −1.3% to +2.7%). The registered values are reproduced exactly (8 of 8) | `P3_06b_cut_bias.tsv` |
| Is the released statistic `net_w` stable when reads are cut to 3 kb? | Cut minus native: −0.016 to +0.004, mean −0.005 (6 lower, 2 higher), all within the registered tolerance. The check is imprecise: 5 draws, and only 91–1,133 sides reach 2 kb after the cut | `P3_06d_netw_cut.tsv` |
| Do lambda control reads sit in the *E. coli* stratum of public runs? | Yes, in 7 of 64 runs (0.04–4.4% of reads; two are counting runs), by a read-length screen. Without them the ARG–MGE net linked fraction changes by at most 0.002 and no p value crosses 0.05. The 8 deep soils show none | `P3_06c_control_strand.tsv`, `P3_06c_control_effect.tsv` |

No registered verdict or estimate is changed by these checks.

One more observation from the measurements: on real runs `mobilink chimeras` flags a median of 12% of the reads it
can judge (4–39%; 1.6% of all reads), against 0.6–2.2% of clean reads in the mock communities. Strains that differ by
an insertion also produce overlap breakpoints. The guide and the draft therefore present the filter as a sensitivity
analysis.

## Release files (step 3.6)

In `pkg/` (local git; commits `a6b44c4`, `ba05843` and the benchmark commit):

| Item | Where |
|---|---|
| Rename to mobilink, version 0.5.0.dev0 | `src/mobilink/`, `pyproject.toml` |
| README, worked example, real-data guide, output columns, method note | `README.md`, `docs/` |
| Example input | `examples/read_lengths.tsv.gz` |
| Helper tables and their build scripts | `resources/` |
| CI workflow (tests on Python 3.10, 3.12, 3.13 and macOS; lint; the worked example; build check) | `.github/workflows/ci.yml` |
| Bioconda recipe draft (three fields to fill at publication) | `recipes/bioconda/meta.yaml` |
| Benchmark: registered texts, frozen scripts, result tables, reports; checksum lists verify | `benchmark/` (built by `scripts/P3_06_benchmark_layout.sh`) |
| Changelog | `CHANGELOG.md` |

- **Name check (2026-11-18):** `mobilink` is free on PyPI, Bioconda and conda-forge
  (`tables/P3/P3_06_name_check.tsv`). The GitHub account `mobilink` exists, and the word is a telecom brand and the
  name of a database-synchronisation product; no bioinformatics tool of that name was found.
- **Old import name.** The analysis scripts of this project import `linkstat` and are unchanged. A project-local
  alias (`compat/linkstat/`, on the path of `envs/linkstat_dev`) maps it to mobilink.
- **Worked example.** Simulated community with a true net linked fraction of 0.125: estimate 0.132 (0.103 to 0.161).
  Reads cut from N50 18.5 kb to 9 and 4 kb: naive share 0.42, 0.31, 0.20; net linked fraction 0.131, 0.130, 0.135.

## Open before publication (yours)

- Repository location and the push; the PyPI and Bioconda releases; the version number at release (now 0.5.0.dev0).
- Licences of the two shipped helper tables (derived from mobileOG-db and BacMet2).
- The frozen benchmark scripts contain the absolute server path.
- Three read-length profiles of gate G3 (30 MB) are not in the repository; they belong in a data archive.
