# Amendment 05: the confirmatory analysis of public data in Phase 3

**Status: approved by the user on 2026-11-16 ("approve both", D56) and frozen the same day, before any confirmatory
run was fetched.** The text is draft v1 with this status block, the approval marks, the file names of the final
sample and the feasibility section; nothing else was changed. The sample was drawn from archive metadata only.
Frozen with: this file's sha256 (`docs/preregistration_amendment_05.sha256`), package commit `9819182` (v0.3.0.dev0),
and the checksummed scripts and inputs in `frozen/amendment05/`.

Scope: hypotheses H1, H2 and H4 of the pre-registration, on public runs that were not used in P0–P2. H3 is reported
descriptively; its between-environment test stays in Phase 4 (amendment 04, 2.6). Our own soils are not public data
and are analysed with the same rules in step 3.4, outside this registration.

## Part A. Sample **[approved 2026-11-16]**

- **Source:** `public/screen/ena_runs.tsv` (ENA portal, fetched 2026-11-14; sha256 `3949b83b…`). Script:
  `scripts/P3_01_select.py`.
- **Eligible run:** ONT; library source METAGENOMIC; strategy WGS (OTHER for PRJNA1414918, as in amendment 02 A2.1);
  selection RANDOM or unspecified; ≥ 1 Gbp; mean read length ≥ 1,500 bp; not a mock, reference or synthetic sample;
  not used in P0–P2, and not from a unit used there.
- **Classes and quotas (48 runs, 12 per class):**

| Class | Project | Eligible runs | Units | Quota |
|---|---|---|---|---|
| Wastewater | PRJNA1207805 (Hong Kong sewage) | 139 | 16 sites | 12 |
| Human gut | PRJNA929328 (transplant patients) | 30 | 23 patients | 6 |
| Human gut | PRJEB108131 (CARMA trial) | 40 | 15 patients | 6 |
| Livestock | PRJNA1279599 (pig manure) | 19 | 19 samples | 4 |
| Livestock | PRJNA1414918 (pig farms) | 15 | 15 farms | 4 |
| Livestock | ERP183232 (chicken caecum) | 32 | 32 birds | 4 |
| Soil | PRJEB58634 (Microflora Danica) | 139 | 138 samples | 4 |
| Soil | PRJEB86780 (UK soils) | 71 | 55 plots | 4 |
| Soil | PRJNA1400213 (Australian soils) | 79 | 79 samples | 4 |

- **Projects of the plan's list that give no eligible run:**
  - PRJNA1088077 (PCR-amplified libraries);
  - ERP188941 (culture-enriched "quasimetagenomes");
  - PRJNA1380421 (selection "other");
  - PRJNA1403057 (mean read length below 1 kb);
  - SRP723315, PRJEB60243 and PRJNA1253088 (no run ≥ 1 Gbp);
  - PRJNA999488 (one run, mean length 1.6 kb, selection "other").

  Hospital wastewater and mine water are therefore not covered.
- **One run per independent unit:** wastewater site, gut patient, pig farm, manure sample, bird, soil plot or sample.
  The run with the most bases represents its unit.
- **Draw:** within a project with quota k, the units are sorted by mean read length and cut into k groups of equal
  size; one unit is drawn at random from each group (seed: CRC32 of the project accession). The sample then spans
  each project's read lengths. The other units, in a seeded random order, are the replacement list.
- **The drawn sample** is `public/screen/p3_selected.tsv` (mean read length 1.7–7.6 kb; identical to the draft's
  `p3_draft.selected.tsv`), with the replacement lists in `public/screen/p3_replacements.tsv`.
- **Replacement** (next unit on the project's list), each recorded: three failed fetches; subsample below 1 Gbp;
  subsample N50 below 2 kb (the registered rule); kit and chemistry both undeterminable.
- **"Kit and basecaller known" (registered rule), read as [approved 2026-11-16]:** the kit is taken from the archive record, the linked
  publication, or the read-end probe (`scripts/P2_05a_kit_probe.py`); the chemistry (R9 or R10) from the record or
  publication. The basecaller is recorded where stated but not required: almost no archive record gives it. R9 runs
  are analysed and reported as their own stratum.
- **Subsample:** 3 Gbp per run as random blocks (48 blocks; 8 larger blocks for runs held only at ENA), or the whole
  run if it has less than 3.3 Gbp. Seeds from the accession, as in P0.

## Part B. Per run

- Annotation and classes as registered: DIAMOND as in the pre-registration; ARG and MRG at identity ≥ 80%; MGE = the
  mobility core at ≥ 60%; subject cover ≥ 80%; same-ORF rule.
- Host strata: Kraken2 on all reads, as in amendment 04 Part 2.1. Unit loci: carrier overlaps, as in Part 2.4.
- Test: the registered net linked test (`linkstat.pairs.net_test_pairs`, package commit `9819182`) for ARG:MGE, ARG:MRG and
  MRG:MGE, with Benjamini–Hochberg across the three pairs within a run and test.
- No chimera filter in the primary analysis (D53).

## Part C. Pooling runs **[approved 2026-11-16]**

- **Pooled test per class and pair:** the runs of a class are analysed as one dataset with the registered test. Host
  strata are run × Kraken2 stratum, so anchors never leave their run and host. Unit loci come from each run's own
  overlaps. Intervals: delete-one-stratum jackknife, as registered.
- Why this and not a meta-analysis of per-run estimates: it uses the registered test unchanged, its p value is exact
  under the clump model, and it needs no per-run standard error (which is unstable in runs with 1–3 effective hosts).
- Limit: the same strain in two units of a class counts as two independent loci. Sensitivity: leave-one-run-out.
- Reported next to it, as descriptive: per-run estimates, and a random-effects summary with its heterogeneity.
- Benjamini–Hochberg across the three pairs within a class and test.

## Part D. Hypotheses and decision rules **[approved 2026-11-16]**

**H1 (the naive share rises with read length).** Naive share = reads with A and B / reads with A, ARG:MGE.
- **H1a, within run (the causal part).** The G1b rule of amendment 02 A3, unchanged, on the new runs: a run counts if
  its subsample N50 is ≥ 9 kb (with the A3 side and interval conditions); its reads are cut to N50 3 kb (20 draws).
  - Supported if at least two-thirds of the counting runs have naive share native / cut ≥ 1.5.
  - Fewer than 6 counting runs: not determined.
- **H1b, across runs (observational).** Spearman correlation between naive share and subsample N50 within each class,
  summed over classes with weights n − 1. One-sided permutation test (N50 permuted within class, 9,999 permutations).
  Supported at p ≤ 0.05. It is reported with the caveat that classes and projects differ in true linkage.

**H2 (the corrected measures do not depend on read length).**
- **Under linkage:** the A3 stability rule on the same counting runs (|native − 3 kb| ≤ max(0.02, 0.25 × |native|),
  intervals overlap), with the frozen Part B estimator. Supported if at least two-thirds of the counting runs meet it.
- **Under no linkage:** the calibration at three read-length profiles is gate G3a (passed). On these runs, the
  Spearman correlation of the ARG:MRG same-molecule estimate with N50 is reported, without a test.

**H4 (ARG–MRG at or below chance in most soils; two-sided).** On the 12 soil runs, same-molecule test, ARG:MRG.
- Per run, with Benjamini–Hochberg across the 12 runs:
  - above chance: q ≤ 0.05 and net_mol > 0;
  - below chance: q ≤ 0.05 and net_mol < 0;
  - at chance: the 90% host-jackknife interval of log(M_A / M_R) lies within ±log 1.25 (the registered TOST);
  - otherwise not determined.
- Pooled over the 12 runs (Part C), the same four outcomes.
- **Supported** if the pooled result is at or below chance and at most 3 of the 12 runs are above chance.
- **Contradicted** if the pooled result is above chance.
- Otherwise **not determined**, and the paper says so with the power available (`linkstat power`).
- "Polluted soils included": the public soils are not polluted. Our own soils (step 3.4) are reported with the same
  rule, outside this registration.

**H3 (descriptive in Phase 3).** Pooled ARG:MGE net_w and net_mol per class, with intervals. No between-class test.

## Part E. Sensitivity analyses (fixed now; none changes a decision above)

- Chimera-flagged reads removed (`linkstat chimeras`, defaults).
- Genus-level host strata; unit distance d = 500 and 2,000 bp.
- ARG and MRG at 60%; MGE at 80%; MGE as the D3 classes; ARG without the intrinsic templates.
- Pooled results without each run in turn, and without runs flagged "few effective hosts".
- R9 runs apart.

## Part F. Consequences

- Each hypothesis is reported as supported, contradicted or not determined, by the rules above. No rule is changed
  after the freeze; any deviation is logged.
- Phase 3 has no stop gate on these results. G4 (submit tier A or continue) is the user's decision.

## Feasibility (development data: the 8 G1b runs; no confirmatory run)

- **Implementation.** Scripts `P3_02_public_run.sh` (one run), `P3_02b_check.py` (status and replacements) and
  `P3_03_confirm.py` (per run, the frozen Part B driver for the A3 columns, pooled test, decision rules). Package
  v0.3.0.dev0 adds M_A / M_R with its 90% host-jackknife interval and the pooling of runs; two G3a datasets reproduce
  exactly, so registered values are unchanged.
- **Dry run of the whole chain** on the 8 G1b runs (3 gut, 3 pig, 2 wastewater), about 25 min on 8 processes:
  - the per-run results equal those of P2 step 2.6 (all p values identical);
  - the A3 columns reproduce G1b (8 of 8 runs count and meet the rule);
  - the run script passed a smoke test on a mock run.
- **Pooled test,** time and memory: pig 3 runs (2.3 M reads) 19 s and 2 GB; gut 3 runs (2.6 M reads) 40 s and 8 GB;
  wastewater 2 runs (1.1 M reads) 211 s and 13 GB. The pooled test is run without the descriptive O/E columns, which
  are slow and memory-hungry with thousands of strata; the per-run tables keep them.
- **What pooling gives** (development data, ARG:MGE net_w): pig 0.067 (0.039 to 0.094) with 92 effective hosts; gut
  0.105 (0.013 to 0.198) with 9; wastewater 0.027 (−0.45 to 0.50) with 1.6, because one run with an E. coli bloom
  carries most units. Such a class shows the flag "few effective hosts"; the leave-one-run-out analysis of Part E
  covers it.
- **H4 machinery:** the 90% intervals of M_A / M_R for ARG:MRG in the development runs were 0.76–0.91, 0.84–1.34 and
  0.30–4.18. "At chance" (within 0.8–1.25) therefore needs many units, and "not determined" is a likely outcome for
  single soil runs. The rule stays as approved.
- **H1a and H2 need runs with N50 ≥ 9 kb.** The sample's mean read lengths are 1.7–7.6 kb, so the counting runs will
  come mostly from wastewater and a few others. Fewer than 6 counting runs gives "not determined", as written.
- **Data pass:** 5 to 35 min of fetching per run, three runs at a time; the host link allows about 7 MB/s in all.
