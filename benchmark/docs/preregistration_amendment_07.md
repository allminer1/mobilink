# Amendment 07: the between-environment analysis of public data in Phase 4

**Status: approved by the user on 2026-11-18 (the reply "approve continue" to the draft and to the proposed gate G5;
D101) and frozen the same day, before any read of the Phase 4 sample was fetched.** The text is draft v1 with this
status block and the approval marks; nothing else was changed. Only run metadata of the archive and the study and run
records of the candidate projects had been read. Frozen with: this file's sha256
(`docs/preregistration_amendment_07.sha256`), package commit `86e3696`, the frozen files of amendments 02 and 05 that
are reused, and the checksummed scripts and inputs in `frozen/amendment07/`. The feasibility section uses development
data (the 64 public runs of Phase 3) and no run of the Phase 4 sample.

Scope: hypothesis H3 of the pre-registration ("ARG–MGE linkage is enriched in wastewater and clinical samples more
than in soil"), whose between-environment test amendment 04 (2.6) left to Phase 4, and three secondary hypotheses.
Every hypothesis below was seen on the Phase 3 runs, so Phase 3 is development data here and the sample is new:
no project with a run fetched in Phases 0–3 is used.

## Part A. Sample **[approved 2026-11-18]**

- **Source:** ENA portal search of 2026-11-18: ONT, library source METAGENOMIC, ≥ 1 Gbp
  (`public/screen/ena_p4_runs_2026-11-18.tsv`, 14,830 runs, sha256 `5dcf7fdf…`).
- **Eligible run:** strategy WGS; selection RANDOM or unspecified; archive mean read length ≥ 2,000; its project has
  no run fetched in Phases 0–3 (28 projects excluded). Script `scripts/P4_01_candidates.py`.
- **Project table** (`public/screen/p4_projects.tsv`, script `scripts/P4_01c_project_table.py`): the archive record
  of every candidate project (217) was read; 100 are included. Excluded: enrichment cultures, laboratory reactors,
  in-vitro and microcosm experiments, mock communities and isolates, amplified libraries, virus-enriched samples,
  adaptive sampling, laboratory and wild animals, samples whose class is not defined or not clear from the record,
  and second projects on a site that is already included.
- **Classes:** human gut (stool); livestock (faeces, rumen content and manure of cattle, pigs, chickens, goats and
  buffalo); wastewater (untreated sewage and influent; effluents left out); activated sludge (full-scale treatment
  plants and digesters); soil; freshwater (lakes, rivers, ponds, groundwater); marine (sea water); sediment.
- **Unit:** the BioSample, or the person, animal, site or plant where the record names it (rules in the project
  table). One run per unit: the run with the longest archive mean read length.
- **Selection:** at most 5 units per project and class, spread over read length by rank, with no random element.
  A class is kept if at least 6 projects and 20 units remain; sediment (9 projects, 18 units) is dropped.
  Script `scripts/P4_01d_units.py`.

| Class | Projects | Runs |
|---|---|---|
| Human gut | 22 | 75 |
| Soil | 25 | 71 |
| Livestock | 10 | 36 |
| Freshwater | 16 | 34 |
| Marine | 12 | 24 |
| Activated sludge | 12 | 21 |
| Wastewater | 10 | 21 |
| **All** | **96** | **282** |

The sample is `public/screen/p4_selected.tsv` (sha256 `0e5f2ff5…`).

- **After the fetch,** a run is rejected if the fetch fails, if its subsample has less than 1 Gbp or a read N50
  below 2 kb, or by the kit rule of amendment 05 (kit and chemistry both unknown from the reads, the record and the
  publication; amplified libraries are not eligible). A rejected run is replaced by the unselected unit of the same
  project and class whose rank by read length is nearest, if there is one. No run is added from another project.
- **A class is analysed** if at least 6 of its projects have at least one accepted run. Otherwise it is reported
  descriptively and takes no part in a test.

## Part B. Per run

- Fetch and annotation: the frozen run script of amendment 05 (`frozen/amendment05/P3_02_public_run.sh`), unchanged:
  3 Gbp as random blocks (the whole run below 3.3 Gbp), DIAMOND against the three databases, Kraken2 on all reads,
  overlaps among the gene-carrying reads, adapter scan. Reads are deleted after annotation.
- Classes of genes: the registered ones (ARG and MRG at ≥ 80% identity, MGE = mobility core at ≥ 60%; subject cover
  ≥ 80%).
- Per-run analysis: the registered net linked test (amendment 04) as in amendment 05 Part B, with the package at the
  commit named in `frozen/amendment07/PACKAGE_COMMIT`, and the Part B cut columns (naive share, N50).

## Part C. Class estimates and their uncertainty **[approved 2026-11-18]**

- **Class estimate.** All accepted runs of a class are pooled with `pool_runs` (the run is part of the stratum, so
  anchors never leave their run and host). The estimate is the pooled ARG:MGE net linked fraction within 2 kb
  (net_w). This is the estimator of amendment 05 Part C.
- **Uncertainty: delete-one-project jackknife.** The pooled estimate is recomputed without each project of the class
  in turn. The standard error is the jackknife standard error over projects, and the interval uses the t
  distribution with (projects − 1) degrees of freedom. Projects, not runs or hosts, are the independent units.
- The jackknife replicates need only the estimate. They are computed with the same package functions in the same
  order as the registered test, without its tests and inner intervals; the script checks for every class that this
  path gives the value of the amendment 05 pooling.
- Reported with each class estimate: the effective number of projects by ARG units, (Σu)² / Σu², where u is the
  number of ARG units of a project. A class in which one project holds most of the units has few effective projects,
  and its interval must be read with that in mind.
- Reported next to it (secondary): the random-effects mean of the project-level estimates (runs of a project pooled;
  DerSimonian–Laird with the Hartung–Knapp interval), for classes in which at least 6 projects have an estimate.

## Part D. Hypotheses and decision rules **[approved 2026-11-18]**

**H3a (primary): ARG–MGE linkage is higher in wastewater and in human gut than in soil.**
- Two contrasts of class estimates (Part C): wastewater − soil, and human gut − soil. "Clinical samples" of the
  registration is read as human gut, the human-associated class of this sample.
- Test: the difference divided by the square root of the sum of the two squared jackknife standard errors, against
  the t distribution with Welch–Satterthwaite degrees of freedom; one-sided; Holm over the two contrasts at 0.05.
- **Supported** if both contrasts are significant. **Partly supported** if one is. **Not supported** if none is.
  **Not determined** if soil, or both other classes, are not analysed (Part A).

**H3b (secondary): the classes differ, and read length does not explain it.**
- Omnibus: the weighted sum of squared deviations of the class estimates from their weighted mean (weights 1 / SE²),
  against the chi-square distribution with (classes − 1) degrees of freedom. Supported at p ≤ 0.05.
- For each pair of analysed classes the contrast is reported with its Holm-adjusted two-sided p value.
- Read length: the same class estimates are computed after the reads of every run are cut to an N50 of 3 kb (one
  breakpoint draw per run, the breakpoint rule of amendment 02; runs below 3 kb are used as they are; each piece
  keeps the host stratum of its read). Reported: the class estimates native and cut, and their rank order. No
  decision depends on it.

**H1c (secondary): within projects, the naive share rises with read length and the net linked fraction does not.**
- Runs of human gut, livestock, wastewater and activated sludge with at least 50 ARG-carrying reads, in projects
  with at least 2 such runs.
- Model: ordinary least squares of the outcome on log2(subsample N50) with one intercept per project and class;
  standard errors clustered by project.
- (i) Naive share (reads with ARG and MGE / reads with ARG): the slope is positive, one-sided p ≤ 0.05.
- (ii) Per-run net_w: the 90% interval of the slope lies within ±0.04 per doubling of N50.
- **Supported** if (i) and (ii) hold. **Not supported** if (i) fails. If (i) holds and the interval of (ii) is not
  within the bounds: **contradicted** if it lies wholly outside them, otherwise **not determined for the net linked
  fraction**.

**H4b (secondary): ARG and MRG on the same read, per class.** The four outcomes of amendment 05 (above chance, below
chance, at chance by the ±log 1.25 bounds, not determined) for the pooled same-molecule result of each analysed
class, with Benjamini–Hochberg across the classes. Reported; no further rule.

**Descriptive, fixed now:** the naive share per class at native length and after the 3 kb cut, next to net_w, and the
rank order of the classes by each; MRG:MGE per class; per-project estimates; the lambda control-read screen of D92 per
run.

## Part E. Sensitivity analyses (fixed now; none changes a decision above)

- Chimera-flagged reads removed (`mobilink chimeras`, defaults).
- Suspected control-strand reads removed (rule of D92).
- Genus-level host strata.
- ARG and MRG at 60%; MGE at 80%; MGE as the broad D3 classes.
- Runs made public before 2021 left out (old chemistry).
- The four projects of one Danish consortium left out together.
- Class estimates with every project given the same weight: the mean of the project-level estimates, for classes in
  which every project has an estimate.

## Part F. Consequences and gate G5 **[approved 2026-11-18]**

- Each hypothesis is reported as supported, partly supported, not supported, contradicted or not determined, by the
  rules above. No rule is changed after the freeze; any deviation is logged.
- The analysis script is run once, after the data pass has ended and the kit rule has been applied.
- **Gate G5 (the user's decision).** Proposed reading of "a clear cross-environment pattern" for a tier B
  submission: H3a is supported; at least 5 classes are analysed; and in the delete-one-project jackknife no single
  project changes the sign of either H3a contrast. If G5 is not met, the fallback is the tier A paper with the
  Phase 4 results as they are.

## Part G. Analyses of Phase 4 outside this registration

These are reported as exploratory, whatever they show: gene classes split further with the existing annotation (MGE by
mobileOG-db category, ARG by drug class); the reanalysis of the 40 datasets of Go and Shrestha 2026, several of which
are in our Phase 3 or Phase 4 samples; projects that sequenced one sample with several extraction or library kits;
and the comparison with isolate or Hi-C data.

## Feasibility (development data: the 64 public runs of Phase 3)

The analysis script (`scripts/P4_03_confirm.py`) was run on the 64 runs of Phase 3 (23 projects; class rule lowered
to 3 projects for this run only, because the development soil class has 3). Output: `tables/P4/dev/confirm/`. These
numbers are not results: the hypotheses were formed on these runs.

| Class | Projects (effective) | Runs | Class estimate (95%, jackknife over projects) | Cut to 3 kb | Naive share | Median N50 |
|---|---|---|---|---|---|---|
| Human gut | 7 (2.8) | 18 | 0.164 (0.071 to 0.257) | 0.172 | 0.47 | 8.8 kb |
| Livestock | 6 (3.1) | 16 | 0.085 (0.037 to 0.133) | 0.098 | 0.25 | 6.0 kb |
| Wastewater | 7 (1.2) | 18 | 0.081 (0.075 to 0.087) | 0.079 | 0.25 | 8.1 kb |
| Soil | 3 (2.7) | 12 | −0.001 (−0.007 to 0.005) | −0.002 | 0.01 | 5.5 kb |

- **The script runs end to end,** and for every class the estimate-only path of the jackknife gives exactly the
  value of the amendment 05 pooling (largest difference 0).
- **H3a on these runs:** wastewater − soil 0.082 (t = 28.1), human gut − soil 0.165 (t = 4.3, Holm p = 0.0025).
  Without any one project, the smallest differences are 0.080 and 0.125. With 6 to 7 projects per class the test has
  room, if the new sample behaves like the old one.
- **One project can dominate a class.** In the development wastewater class one project holds 93% of the ARG units
  (1.2 effective projects), and the interval is narrow because the other six projects happen to agree with it. The
  Phase 4 sample takes at most 5 runs per project, which limits this; the effective number of projects is reported
  with every class, and the random-effects mean of the project-level estimates is reported next to the pooled one
  (here: human gut 0.097, livestock 0.058, wastewater 0.092).
- **The naive share orders the classes differently:** human gut > wastewater > livestock > soil, against human gut >
  livestock > wastewater > soil by the net linked fraction, native and cut to 3 kb. The wastewater runs have longer
  reads than the livestock runs.
- **H1c on these runs:** naive share +0.089 per doubling of N50 (one-sided p = 0.014; 38 runs, 6 projects); net_w
  +0.019 (90% interval −0.008 to 0.046; 36 runs, 5 projects), which is "not determined for the net linked fraction"
  by the rule. The Phase 4 sample has about 60 projects in the four classes of H1c.
- **Sparse classes:** only 1 of 3 development soil projects has a project-level estimate, and 4 of 7 wastewater
  projects (three are single runs of sludge with few ARG). The pooled class estimate exists for every class. This is
  why the pooled estimate with a project jackknife is primary and the mean of project-level estimates is secondary.
- **Cost:** about 1.2 hours for 64 runs on 8 processes under load; the Phase 4 run is estimated at 6 to 10 hours.
  Two functions of the package were made faster for this (package commits `05cdfde` and `86e3696`); three registered
  per-run tables are reproduced exactly with them (58 columns, difference 0) and all 104 tests pass.
