# Benchmark and pre-registered validation

This folder holds the scripts, the registered texts and the result tables behind every validation claim of mobilink.
The files are unchanged copies from the project in which the method was developed. The folder mirrors the layout of
that project, so the checksum lists still verify:

```bash
cd benchmark
for d in amendment03 amendment04 amendment05 amendment06; do (cd frozen/$d && sha256sum -c --quiet SHA256SUMS); done
(cd frozen/partB && sha256sum -c --quiet ../partB.sha256)
sha256sum -c docs/preregistration_amendment_05.sha256 docs/preregistration_amendment_06.sha256
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

Deviations from the registered rules are logged in `docs/preregistration_amendment_05_deviations.md` and
`docs/preregistration_amendment_06_deviations.md`.

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
| Figures and supplementary tables of the paper | `scripts/P3_07_figures.py`, `P3_07_supp_tables.py` | | |

## Package versions

The package was called `linkstat` while these analyses ran. Each frozen folder names the commit it used
(`frozen/*/PACKAGE_COMMIT`); all of them are in the history of this repository.

| Analysis | Commit | Version |
|---|---|---|
| G1b | none: the estimator is the frozen script `frozen/partB/P1_03_estimators.py` | |
| G2 | `7b665ca` | 0.1 dev |
| G3a, G3b | `100c199` | 0.2 dev |
| Stage 1, stage 2 | `9819182` | 0.3.0.dev0 |

The release (0.5) gives the same numbers as commit `9819182` for the registered per-run analysis; see
`tables/P3/P3_06_release_check.tsv`.

## Running the scripts again

1. Check out the commit of the analysis and install it (`pip install -e .`). For the release instead, install
   mobilink and put `benchmark/compat` on `PYTHONPATH`: it makes `import linkstat` resolve to mobilink.
2. The scripts expect the project root in the variable `L` at the top of each file, with the folders `public/`,
   `db/`, `tables/` and `tmp/` under it. Set it to your own path.
3. The public reads are fetched from SRA or ENA by the run scripts and deleted after annotation. The accession lists
   are in `public/screen/` (`p3_final.tsv` for stage 1, `p3x_final.tsv` for stage 2). Stage 2 (25 runs tried, 3 Gbp
   each) took about 12 hours with 36 threads on a busy server.
4. A registered analysis is run once. If you repeat it, write to a new output folder; the drivers refuse to
   overwrite an existing one.

The simulations (G2, G3) need no download except the three read-length profiles named above.
