# Benchmark and pre-registered validation

This folder holds the scripts, the registered texts and the result tables behind every validation claim of mobilink.
The files are unchanged copies from the project in which the method was developed. The folder mirrors the layout of
that project, so the checksum lists still verify:

```bash
cd benchmark
for d in amendment03 amendment04 amendment05 amendment06 amendment07; do (cd frozen/$d && sha256sum -c --quiet SHA256SUMS); done
(cd frozen/partB && sha256sum -c --quiet ../partB.sha256)
sha256sum -c docs/preregistration_amendment_05.sha256 docs/preregistration_amendment_06.sha256 docs/preregistration_amendment_07.sha256
```

Three files named in `frozen/amendment04/SHA256SUMS` are not in the repository because of their size (read-length
profiles, 30 MB together): `lengths/279.lengths.tsv.gz`, `public/SRR32022521/SRR32022521.lengths.tsv.gz` and
`public/SRR23278673/SRR23278673.lengths.tsv.gz`. They belong in the data archive of the paper.

## What was tested, and where the result is

Every gate and hypothesis was fixed in writing before its data were analysed (`docs/`). Failed gates are kept.

| Step | Question | Result | Registered in | Scripts | Tables | Report |
|---|---|---|---|---|---|---|
| G1b | Does the naive share change with read length while the net linked fraction stays stable, in public runs? | **Pass**, 8 of 8 runs in 3 environments | amendment 02 | `frozen/partB/`, `frozen/partB_gate_P1_08_g1b_gate.py`, `scripts/P1_07_g1b_driver.sh` | `tables/P1_07_g1b_MC.tsv`, `tables/P1_08_g1b_gate.tsv` | `reports/P1/G1b_memo.md` |
| G2 | Is the analytic read-level O/E test calibrated? | **Fail** (G2a: 4 of 7 cases within 5%; G2b: 3 of 6 simulation settings pass). The net linked test replaced it as the read-level test | amendment 03 | `frozen/amendment03/`, `scripts/P1_15_g2_driver.sh` | `tables/G2/` | `reports/P1/G2_memo.md` |
| G3a | Is the net linked test calibrated, and do the intervals cover the truth? | **Pass**: size 5.0% and 4.4% at a nominal 5% (9,000 datasets without linkage); lowest coverage 0.914 | amendment 04 | `frozen/amendment04/`, `scripts/G3a_run.sh` | `tables/G3/G3a_*` | `reports/P2/G3_report.md` |
| G3b | Is the false discovery rate controlled when several pairs are tested (Benjamini–Hochberg)? | **Pass**: mean false share 0.031 (`net_w`) and 0.029 (same molecule) | amendment 04 | `frozen/amendment04/G3b_summary.py` | `tables/G3/G3b_*` | `reports/P2/G3_report.md` |
| Stage 1 | H1a, H2: cut reads to 3 kb in public runs. H1b, H4 | H1a and H2 **not determined** (3 of 3 counting runs met the criterion, 6 were needed); H1b supported (p = 0.0475); H4 not determined | amendment 05 | `frozen/amendment05/`, `scripts/P3_02_driver.sh`, `scripts/P3_03_run.sh` | `tables/P3/confirm/` | `reports/P3/P3_03_confirmatory.md` |
| Stage 2 | H1a, H2 again, on runs chosen for a read N50 of at least 9 kb after stage 1 was known | H1a and H2 **supported**, 11 of 11 counting runs | amendment 06 | `frozen/amendment06/`, `scripts/P3X_*` | `tables/P3/confirm_ext/` | `reports/P3/P3_03x_stage2.md` |
| Survey | H3a: is ARG–MGE linkage stronger in wastewater and in human gut than in soil? H3b: do seven environment classes differ? H1c: within projects, does the naive share follow read N50 while the net linked fraction does not? H4b: ARG–MRG per class. Gate G5: is the cross-environment pattern clear? | 241 runs, 86 projects. H3a **not supported** (differences −0.020 and +0.016, SE 0.060). H3b supported (p = 0.006). H1c: naive share slope +0.117 per doubling (supported); net linked fraction +0.016, 90% interval −0.011 to +0.044, **not determined** against the bound of ±0.04. H4b above chance in activated sludge only. **Gate G5 not met** | amendment 07 | `frozen/amendment07/`, `scripts/P4_0[1-3]*`, `P4_freeze_check.sh` | `tables/P4/confirm/` | `reports/P4/P4_03_confirmatory.md` |

Deviations from the registered rules, and notes on their application, are logged in
`docs/preregistration_amendment_05_deviations.md`, `docs/preregistration_amendment_06_deviations.md` and
`docs/preregistration_amendment_07_deviations.md`.

Outside the registration (exploratory, reported as such):

| Analysis | Scripts | Tables | Report |
|---|---|---|---|
| Sensitivity of stage 1 (chimeras, class definitions, host rank) | `scripts/P3_04_sensitivity.py` | `tables/P3/P3_04_sensitivity.tsv` | `reports/P3/P3_03_confirmatory.md` |
| Mock communities (chimera flags) | | | `reports/P2/P2_05_mocks_report.md` |
| Eight deep soils, read level | | `tables/P3/P3_04_soil_reads.tsv` | `reports/P3/P3_04_soils_read_level.md` |
| 27 low-coverage soil libraries | `scripts/P3_04d_s1_*` | `tables/P3/P3_04d_s1.*` | `reports/P3/P3_04d_low_coverage.md` |
| Contigs: label permutation against the mobilink test, with a no-linkage simulation | `scripts/P3_05a_pilot_loci.py`, `P3_05b_pilot_test.py`, `P3_05d_soil_contig_test.py`, `P3_05e_null_sim.py` | `tables/P3/P3_05_*` | `reports/P3/P3_05_contigs.md` |
| Run time and memory of the released commands; release check against the registered tables (exact in 64 of 64 runs) | `scripts/P3_06_runtime*.sh`, `P3_06_runtime_table.py` | `tables/P3/P3_06_runtime*.tsv`, `P3_06_release_check.tsv` | |
| Why the cut estimate is slightly lower; the unit-level estimator on cut reads | `scripts/P3_06b_cut_bias.py`, `P3_06d_netw_cut.py` | `tables/P3/P3_06b_cut_bias.tsv`, `P3_06d_netw_cut.tsv` | |
| Lambda control reads in the public runs and their effect | `scripts/P3_06c_control_*` | `tables/P3/P3_06c_control_*.tsv` | |
| Figures and supplementary tables of the paper | `scripts/P3_07_figures.py`, `P3_07_supp_tables.py`, `P4_08_figures.py`, `P4_09_supp_tables.py` | | |
| Survey: sensitivity analyses fixed with amendment 07 (samples, chimeras, control reads, genus strata, class definitions) | `scripts/P4_04_sens_samples.py`, `P4_04b_sens_variants.py` | `tables/P4/sens/` | `reports/P4/P4_03_confirmatory.md` |
| Survey: sub-classes of ARG (mechanism, drug family) and of MGE (integration/excision, transfer) | `scripts/P4_06a_*`, `P4_06b_*`, `P4_06c_*` | `tables/P4/sub/` | `reports/P4/P4_06_subclasses.md` |
| A published comparison of eight environments (Go and Shrestha 2026): their calls against read length, and the same 40 runs with this pipeline | `scripts/P4_07_driver.sh`, `P4_07a_their_calls.py`, `P4_07c_estimates.py`, `P4_07e_compare.py` | `tables/P4/reanalysis/` | `reports/P4/P4_07_reanalysis.md` |
| One sample sequenced after several extraction or library methods (three public sets, 88 runs) | `scripts/P4_07b_driver.sh`, `P4_07c_estimates.py`, `P4_07d_kit_slopes.py` | `tables/P4/kits/` | `reports/P4/P4_07_reanalysis.md` |
| Version 1.0: release check against the registered per-run tables (exact in 64 of 64 runs), run time and memory, and `mobilink meta` against the registered class estimates (exact) | `scripts/P4_09_release*.sh`, `P4_09_release_table.py` | `tables/P4/v1check/` | `reports/P4/P4_09_package_v1.md` |

## Package versions

The package was called `linkstat` while these analyses ran. Each frozen folder names the commit it used
(`frozen/*/PACKAGE_COMMIT`).

**Commit numbers.** On 2026-11-23, before the repository was made public, the author field and the messages of its
commits were edited. No file changed, but every commit got a new number. The registered texts, the `PACKAGE_COMMIT` files and the
reports are checksummed or kept as written, so they still give the numbers from before. `COMMITS.tsv` lists the old
and the new number of every commit, and each old number that is cited has a tag `orig-<old number>` on the commit
with the same content, for example `git checkout orig-86e3696`. The freeze-check scripts compare the checked-out
commit with `PACKAGE_COMMIT` and therefore stop on the new numbers; they belong to analyses that were run once and
are not meant to be run again.

| Analysis | Commit as cited in the frozen files | Version |
|---|---|---|
| G1b | none: the estimator is the frozen script `frozen/partB/P1_03_estimators.py` | |
| G2 | `7b665ca` (tag `orig-7b665ca`) | 0.1 dev |
| G3a, G3b | `100c199` (tag `orig-100c199`) | 0.2 dev |
| Stage 1, stage 2 | `9819182` (tag `orig-9819182`) | 0.3.0.dev0 |
| Survey (amendment 07) and all analyses of Phase 4 | `86e3696` (tag `orig-86e3696`) | 0.5.0.dev0 |

Version 0.5 gives the same numbers as commit `9819182` for the registered per-run analysis
(`tables/P3/P3_06_release_check.tsv`), and so does version 1.0 (`tables/P4/v1check/P4_09_release_check.tsv`).
Version 1.0 also reproduces the class estimates, jackknife replicates and contrasts of the survey with
`mobilink meta` (`tables/P4/v1check/meta_ARG_MGE.*.tsv`). Its chimera flags differ from version 0.5 by a correction
of ties (0.3% of the flags removed, none added; see `CHANGELOG.md`); the chimera filter is a sensitivity analysis
only and was run with version 0.5.

## Running the scripts again

1. Check out the commit of the analysis and install it (`pip install -e .`). For the release instead, install
   mobilink and put `benchmark/compat` on `PYTHONPATH`: it makes `import linkstat` resolve to mobilink.
2. The scripts expect the project root in the variable `L` at the top of each file, with the folders `public/`,
   `db/`, `tables/` and `tmp/` under it. Set it to your own path.
3. The public reads are fetched from SRA or ENA by the run scripts and deleted after annotation. The accession lists
   are in `public/screen/` (`p3_final.tsv` for stage 1, `p3x_final.tsv` for stage 2, `p4_final.tsv` for the survey,
   `p4_07_runs.tsv` and `p4_07b_runs.tsv` for the published comparison and the kit sets). Stage 2 (25 runs tried, 3 Gbp
   each) took about 12 hours with 36 threads on a busy server.
4. A registered analysis is run once. If you repeat it, write to a new output folder; the drivers refuse to
   overwrite an existing one.

The simulations (G2, G3) need no download except the three read-length profiles named above.

## What is not in this folder

- Reads, and the per-run tables made from them (gene calls, read lengths, host labels, overlaps): about 100 GB. They
  belong in the data archive of the paper.
- Databases (CARD, BacMet, mobileOG-db, Kraken2). Only the Kraken2 taxonomy table and the list of mobility-core
  entries are here, because checksum lists name them.
- The raw archive records of the 217 candidate projects of the survey; their digest is in
  `public/screen/p4_records.digest.tsv`.
- The read-level table of Go and Shrestha (Zenodo, doi:10.5281/zenodo.23492757) and the copy of its read rows made
  by `scripts/P4_07a_their_calls.py`. Download the table from Zenodo to run that script.

The reports and the registered texts were written during the project and are kept as they were. They name decisions
by number (D1, D2, ...); the decision log is part of the project archive. "The user" in these texts is the
investigator who led the project and took the decisions that the plans reserve for a person (gates, changes to
registered definitions, anything public).
