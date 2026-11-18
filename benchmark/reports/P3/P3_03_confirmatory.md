# P3 step 3.3: confirmatory analysis (amendment 05)

2026-11-17. Frozen script `frozen/amendment05/P3_03_confirm.py`, package commit `9819182`, run once on the final
list of 48 public runs (`tables/P3/confirm/`, verdicts in `verdict.txt`). The freeze check passed before the run.
Rules are those of amendment 05 Part D; none was changed.

## Result

**Two of the four tested hypotheses could not be decided, because the sample has too few long-read runs. One is
weakly supported, and the soil hypothesis cannot be tested with the registered classes.**

| Hypothesis | Registered verdict | Basis |
|---|---|---|
| H1a: naive share falls when the reads of a run are cut | **not determined** (stage 2: supported, 11 of 11) | 3 counting runs (6 needed); all 3 meet the criterion |
| H1b: naive share rises with read N50 across runs | **supported**, weakly | one-sided permutation p = 0.0475 |
| H2: net linked fraction does not change when reads are cut | **not determined** (stage 2: supported, 11 of 11) | the same 3 counting runs; all 3 are stable |
| H4: ARG–MRG at or below chance in soil | **not determined** | no ARG–MRG read in 27.8 Gbp; 0.67 expected by chance |
| H3: ARG–MGE linkage per class | descriptive | net linked fraction within 2 kb: gut 0.21, livestock 0.10, wastewater 0.08, soil 0.00 |

- **H1a and H2 fail on sample size, not on the data.** A run counts only if its N50 is ≥ 9 kb. Four of the 48 runs
  reach that, and one of them (soil) has too few ARG.
- **Every run that can be evaluated points the same way.** Of the 43 runs with both estimates, the net linked
  fraction is stable in 42. Of the 25 with N50 ≥ 6 kb, the naive share falls ≥ 1.5-fold in all 25. This is
  descriptive: the registered rule counts only runs with N50 ≥ 9 kb.
- **Decided 2026-11-17 (D70, D72):** an extension for H1a and H2 is registered as amendment 06 (a second stage of
  runs with N50 ≥ 9 kb, chosen by a fixed rule before any of them was looked at). See "What this leaves open".
- **Sensitivity analyses (section 6):** no conclusion changes. ARG:MGE stays above chance in gut, livestock and
  wastewater in every variant; the low ARG:MRG counts within 2 kb are not robust to the identity threshold.

## 1. Sample

| Class | Runs | Projects | Subsample (Gbp) | N50 (kb), min–median–max | Runs with N50 ≥ 9 kb | Chemistry |
|---|---|---|---|---|---|---|
| Human gut | 12 | 1 | 34.9 | 4.7–7.7–12.5 | 3 | R9.4 |
| Livestock | 12 | 3 | 33.2 | 3.2–5.7–8.3 | 0 | R10.4.1 (4 runs); not stated (8 runs) |
| Soil | 12 | 3 | 27.8 | 4.1–5.5–10.4 | 1 | R10.4.1 |
| Wastewater | 12 | 1 | 36.0 | 5.3–7.8–8.3 | 0 | R10.4.1 |

- 56 runs were fetched; 8 were excluded and replaced (`docs/preregistration_amendment_05_deviations.md`):
  6 gut runs made from PCR-amplified libraries (deviation 1) and 2 chicken runs with no kit or chemistry on record.
- All 48 runs have ≥ 1 Gbp and N50 ≥ 2 kb. No fetch failed.
- The gut class is one project and all R9.4, so chemistry and class cannot be separated there.

## 2. H1: the naive share and read length

**H1a (within run).** Reads of each run were cut to N50 3 kb (20 draws); naive share = reads with ARG and MGE / reads
with ARG.

| N50 of the run | Runs | Both estimates | Naive share native / 3 kb ≥ 1.5 | Median ratio (range) |
|---|---|---|---|---|
| ≥ 9 kb (the registered counting range) | 4 | 3 | 3 | 2.00 (1.81–3.41) |
| 6 to < 9 kb | 25 | 22 | 22 | 1.91 (1.50–8.47) |
| < 6 kb | 19 | 18 | 7 | 1.40 (1.07–3.62) |

- Registered: 3 counting runs, all gut, all ≥ 1.5. Fewer than 6, so **not determined**.
- The ratio grows with the N50 of the run (Spearman 0.81 over the 36 non-soil runs), as expected when the cut is to
  a fixed 3 kb. Runs below 6 kb are cut by less than half and mostly stay under 1.5.

**H1b (across runs).** Within-class Spearman correlation of naive share with N50: wastewater 0.66, gut 0.39,
livestock 0.08, soil −0.13. Weighted sum 11.02, one-sided permutation p = 0.0475: **supported**, at the threshold.
Classes and projects differ in true linkage, so this is the weaker, observational half of H1.

## 3. H2: the net linked fraction and read length

- Registered (under linkage): the 3 counting runs are all stable with overlapping intervals; **not determined**
  for the same reason as H1a.

| Counting run | N50 (kb) | Naive share: native → 3 kb | Net linked: native → 3 kb → 2 kb |
|---|---|---|---|
| SRR23278716 | 12.5 | 0.667 → 0.334 | 0.283 → 0.279 → 0.273 |
| SRR23278658 | 9.0 | 0.513 → 0.284 | 0.239 → 0.234 → 0.237 |
| SRR23278671 | 12.3 | 0.446 → 0.131 | 0.064 → 0.063 → 0.066 |

- All runs (descriptive): 42 of 43 runs with both estimates are stable; the absolute change is 0.0009 at the
  median and 0.0063 at most (1.6% of the estimate at the median). The one failure is a soil run whose 3 kb interval
  could not be computed.
- Under no linkage (descriptive, as registered): Spearman of the ARG:MRG same-molecule estimate with N50 over 26
  runs is 0.27 (p = 0.18).

## 4. H4: ARG–MRG in soil

- The 12 soil runs hold 770 ARG and 540 MRG units in 8.56 million reads. No read carries both; 0.67 are expected
  by chance on the same read and 0.24 within 2 kb.
- Per run: 12 of 12 "not determined" (no events). Pooled: "not determined". Verdict: **not determined**.
- **Power** (`linkstat power`, soil run ERR11523675, the one with the most ARG): at 3 Gbp, 0.08 co-carrying reads
  are expected. A depletion to 0.14 of chance becomes detectable at about 300 Gbp, and to 0.5 at about 1,000 Gbp.
  Read-level data of this kind cannot test H4 with the registered classes.
- The contig level is the place for this question; see `reports/P3/P3_05_contigs.md` (our soils, outside the
  registration: no depletion, all "not determined").

## 5. H3 and the other pairs (descriptive)

**ARG:MGE, pooled per class.**

| Class | ARG units | Net linked within 2 kb (95%) | Net linked, same read (95%) | Share ratio M_A / M_R (90%) | Runs above chance within 2 kb |
|---|---|---|---|---|---|
| Human gut | 12,830 | 0.214 (0.179 to 0.250) | 0.508 (0.442 to 0.574) | 2.93 (2.46–3.50) | 7 of 12 |
| Livestock | 17,090 | 0.102 (0.086 to 0.118) | 0.261 (0.226 to 0.296) | 2.58 (2.35–2.83) | 5 of 12 |
| Wastewater | 22,920 | 0.080 (0.065 to 0.096) | 0.184 (0.144 to 0.225) | 2.20 (1.95–2.48) | 11 of 12 |
| Soil | 770 | −0.001 (−0.004 to 0.002) | 0.001 (−0.014 to 0.015) | 1.05 (0.43–2.55) | 0 of 12 (10 without an estimate) |

- q ≤ 0.001 in gut, livestock and wastewater for both measures; q = 0.97 and 0.95 in soil.
- Runs differ within a class (random-effects summary, `tables/P3/P3_03.random_effects.tsv`): gut 0.17 (0.09 to
  0.25), I² 97%; livestock 0.10 (0.07 to 0.13), I² 82%; wastewater 0.09 (0.07 to 0.11), I² 70%.
- No test between classes is made (as registered).

**ARG:MRG, pooled per class (observed / expected events).**

| Class | Within 2 kb | p | Same read | p | Share ratio (90%) | By the H4 rule |
|---|---|---|---|---|---|---|
| Human gut | 0 / 16.5 | 0.051 | 73 / 81.2 | 0.92 | 1.21 (0.73–1.98) | not determined |
| Livestock | 0 / 13.3 | 0.013 | 46 / 48.0 | 0.98 | 1.04 (0.84–1.28) | not determined |
| Wastewater | 75 / 149.2 | 0.22 | 616 / 508.0 | 0.26 | 1.26 (1.12–1.42) | not determined |
| Soil | 0 / 0.24 | 0.89 | 0 / 0.67 | 0.74 | no estimate | not determined |

- On the same read, ARG and MRG meet about as often as chance predicts in gut and livestock. Livestock misses "at
  chance" narrowly (upper bound 1.28 against 1.25).
- Within 2 kb the counts are below expectation in three classes, but the expectation sits in very few loci
  (0.4, 1.6 and 5.4 effective events; the tool flags all three). In gut and livestock one run carries most of it.
  This is not evidence of depletion.
- **MRG:MGE** is above chance in wastewater (within 2 kb 1,258 / 548, q = 0.0002; net linked 0.027) and not
  significant elsewhere.

## 6. Sensitivity analyses (Part E)

Finished 2026-11-17 15:28 UTC (`scripts/P3_04_sensitivity.py` → `tables/P3/P3_04_sensitivity.tsv`, 93 pooled
analyses). As registered, none of them changes a verdict above. The row "primary" reproduces the frozen analysis.

**ARG:MGE, net linked within 2 kb, pooled per class.**

| Class | Primary (95%) | Other Part E variants, lowest to highest | Chimera-flagged reads removed | MGE as the D3 classes | Without each run in turn |
|---|---|---|---|---|---|
| Human gut | 0.214 (0.179 to 0.250) | 0.181 (MGE at 80%) to 0.221 (genus strata) | 0.207 | 0.352 | 0.192 to 0.236 |
| Livestock | 0.102 (0.086 to 0.118) | 0.078 (ARG at 60%) to 0.106 (genus strata) | 0.108 | 0.307 | 0.093 to 0.110 |
| Wastewater | 0.080 (0.065 to 0.096) | 0.046 (ARG at 60%) to 0.083 (genus strata) | 0.059 | 0.169 | 0.076 to 0.087 |
| Soil | −0.001 (−0.004 to 0.002) | −0.002 to 0.001 | 0.000 | −0.006 | −0.002 to 0.000 |

- Gut, livestock and wastewater stay above chance in every variant and without any single run (q ≤ 0.006 in all).
  Soil stays at chance in every variant (q ≥ 0.23).
- The size depends on the class definitions more than on anything else. The broad D3 MGE classes give about twice
  to three times the primary value; ARG at 60% identity lowers it, most in wastewater (0.046).
- Removing chimera-flagged reads (63,000 to 229,000 reads per class) lowers the wastewater value by about a quarter
  (0.080 to 0.059) and leaves gut and livestock within their intervals.
- One soil interval excludes zero (d = 2,000: −0.0031 to −0.0018), on 2 observed events against 2.3 expected
  (q = 0.97). This is not evidence of depletion.
- Runs flagged "few effective hosts": one livestock run; without it the value is 0.101.
- R9 runs apart: all 12 gut runs are R9.4 and no run of another class is known to be R9, so the gut row is the R9
  stratum and no further split exists.

**ARG:MRG.**

- **Within 2 kb, classes at 80%.** The count stays below expectation in every variant: gut 0 observed against 15.2
  to 17.0 expected, livestock 0 against 11.1 to 13.7, wastewater 64 to 75 against 141 to 152. The effective events
  stay at 0.4 to 0.7 (gut), 1.3 to 1.6 (livestock) and 5.2 to 5.5 (wastewater).
- **Within 2 kb, classes at 60%.** The direction changes: gut 12 against 28.9 (q = 0.45), livestock 44 against 23.4
  (q = 0.15), wastewater 373 against 269.3 (q = 0.20), soil 15 against 5.6 (q = 0.15). The low counts at 80% are
  therefore not robust to the identity threshold.
- **Without each run in turn.** Gut: without SRR23278658 the expectation falls from 16.5 to 5.8 (p = 0.12).
  Livestock: without ERR17169904 it falls from 13.3 to 5.1 (p = 0.044). This confirms that one run carries most of
  the expectation in each class.
- **Same read.** Gut and livestock stay near chance (share ratio 1.02 to 1.65 and 0.80 to 1.16 without each run in
  turn). Wastewater stays modestly above chance: share ratio 1.19 to 1.29 across the 80% variants, all 90% intervals
  above 1, but the clump test is not significant (q = 0.22 to 0.26). With classes at 60% it is significant (1,305
  against 927; ratio 1.38, 1.28 to 1.48; q = 0.004).
- **Soil (H4).** No estimate in any variant at 80% (0 observed against 0.4 to 0.7 expected on the same read). At 60%:
  28 against 16.3 on the same read, share ratio 1.36 (0.92 to 2.02), not determined. H4 stays "not determined" in
  every variant.

**MRG:MGE.**

- Wastewater is above chance in every variant and without any single run (net linked 0.017 to 0.030; 0.068 with the
  D3 classes; q ≤ 0.0004).
- Gut (112 against 60.0, q = 0.051, 4 effective hosts) and livestock (61 against 35.4, q = 0.11) are not robust: q
  runs from 0.006 to 0.18 in gut and from 0.046 to 0.68 in livestock across the variants.
- Soil is at chance at 80%. At 60% the same-read share ratio is 0.73 (0.56 to 0.94), with q = 0.16.

**Outside Part E: the six PCR-amplified gut runs of deviation 1** (not registered; shown for transparency).

- ARG:MGE with the six runs added to the 12: 0.199 (0.161 to 0.237), against 0.214 without them.
- The six runs alone: 0.084 (−0.055 to 0.222) with 5 effective hosts; 1,118 events within 2 kb against 189 expected
  (q = 0.0003). Four runs give an estimate on their own (0.038, 0.052, 0.141, 0.192); two have too few sides.
- ARG:MRG in the six runs alone: 3 against 67.0 within 2 kb (1.3 effective events) and 73 against 108.9 on the same
  read (share ratio 0.65, 0.17 to 2.41). Two runs carry almost all of the expectation.
- The amplified runs do not change any conclusion. They stay out of the registered sample.

## What this leaves open

1. **H1a and H2.** These are the paper's main claim. The earlier gate G1b passed 8 of 8 on development runs, and all
   25 confirmatory runs with N50 ≥ 6 kb agree, but the registered count is 3 of 3 with 6 needed. 
   **Decided 2026-11-17 (D70, D72): an extension is registered as amendment 06.** It adds runs with N50 ≥ 9 kb chosen by a fixed
   rule, keeps the Part D rule unchanged and reports the result as a second stage; the stage-1 verdict stands. The
   stage-2 result goes to `reports/P3/P3_03x_stage2.md`.
   **Stage 2, run once on 2026-11-17/18: H1a supported and H2 supported, each in 11 of 11 counting runs** (16 runs
   analysed: gut 6, livestock 4, wastewater 6). Both stages together: 14 of 14.
2. **H4** cannot be tested on reads with the registered classes at any realistic depth. The paper should say so and
   carry the ARG–MRG question at the contig level.
3. Feasibility note of amendment 05 expected the counting runs to come mostly from wastewater. The wastewater
   subsamples reached N50 5.3–8.3 kb, so none counted.

## Files

- `tables/P3/confirm/`: `verdict.txt`, `verdict.tsv`, `per_run.tsv`, `partB.tsv`, `pooled.tsv`, `H4_runs.tsv`, and
  copies of the run list and status.
- `tables/P3/P3_03.random_effects.tsv`, `tables/P3/P3_03.partB_all_runs.tsv` (`scripts/P3_03b_descriptive.py`).
- `tables/P3/P3_04_sensitivity.tsv` (`scripts/P3_04_sensitivity.py`; log `logs/P3/P3_04.log`).
- `public/screen/p3_kit_probe.tsv`, `p3_project_kits.tsv`, `p3_manual_exclusions.tsv`.
