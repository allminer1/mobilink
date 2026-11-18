# P3 step 3.3x: stage 2 for H1a and H2 (amendment 06)

2026-11-18. Rules: amendment 06 (frozen 2026-11-17), which applies the Part D rule of amendment 05 unchanged to newly
selected runs. Code: the frozen stage-1 script and package commit `9819182`. The analysis was run once
(`scripts/P3X_03_run.sh`, 2026-11-17 23:42 to 2026-11-18 00:03 UTC) and will not be rerun. Tables:
`tables/P3/confirm_ext/` (`verdict.txt`, `partB.tsv`, `per_run.tsv`, `pooled.tsv`).

## Verdicts

| Hypothesis | Stage 2 verdict | Basis |
|---|---|---|
| H1a: naive share falls when the reads of a run are cut | **supported** | 11 of 11 counting runs have a ratio ≥ 1.5 (median 3.61, range 1.99–7.66) |
| H2 (under linkage): net linked fraction does not change when reads are cut | **supported** | 11 of 11 counting runs are stable with overlapping intervals (difference: median 0.0014, largest 0.0063) |

Wording fixed in amendment 06, for the paper:

- **H1a.** "Stage 1 (amendment 05): not determined, 3 of 3 counting runs met the criterion, 6 were needed. Stage 2
  (amendment 06, runs chosen for N50 ≥ 9 kb after stage 1 was known): supported, 11 of 11 counting runs."
- **H2.** "Stage 1 (amendment 05): not determined, 3 of 3 counting runs met the criterion, 6 were needed. Stage 2
  (amendment 06, runs chosen for N50 ≥ 9 kb after stage 1 was known): supported, 11 of 11 counting runs."

The stage-1 verdicts are not replaced. H1b, H3 and H4 were not tested again; their stage-1 results stand.

## Counting runs

A run counts if its subsample N50 is ≥ 9 kb, at least 30 ARG sides reach 2 kb in the native and the cut version, and the
half-width of the interval at 3 kb is ≤ 0.10. Pair ARG:MGE, window 2 kb; "3 kb" is the same reads cut to an N50 of 3 kb
(20 draws).

| Run | Class | Kit; chemistry | N50 | Naive share: native → 3 kb (ratio) | Net linked: native (interval) | Net linked: 3 kb (interval) | Difference |
|---|---|---|---|---|---|---|---|
| SRR29980985 | human gut | ligation (record); R10.4.1 | 18,536 | 0.472 → 0.151 (3.13) | 0.085 (0.074 to 0.097) | 0.080 (0.059 to 0.101) | +0.0044 |
| SRR18491323 | human gut | ligation, native barcodes (record); R9.4 | 13,201 | 0.393 → 0.083 (4.74) | 0.017 (0.010 to 0.025) | 0.016 (-0.001 to 0.032) | +0.0010 |
| SRR33225186 | human gut | native barcoding (record); R10.4.1 | 11,990 | 0.541 → 0.248 (2.18) | 0.153 (0.146 to 0.160) | 0.151 (0.135 to 0.167) | +0.0022 |
| DRR351862 | human gut | native barcoding (reads); not stated | 13,477 | 0.450 → 0.210 (2.14) | 0.111 (0.091 to 0.131) | 0.105 (0.068 to 0.141) | +0.0063 |
| ERR15826294 | human gut | ligation (reads); not stated | 10,771 | 0.171 → 0.047 (3.61) | 0.024 (0.014 to 0.036) | 0.023 (-0.000 to 0.045) | +0.0014 |
| SRR16329124 | human gut | ligation (reads); not stated | 14,904 | 0.334 → 0.167 (1.99) | 0.143 (0.130 to 0.157) | 0.138 (0.109 to 0.167) | +0.0047 |
| ERR16013023 | livestock | native barcoding (reads); not stated | 11,136 | 0.035 → 0.005 (7.66) | -0.002 (-0.002 to -0.001) | -0.002 (-0.003 to -0.001) | -0.0001 |
| SRR25823913 | livestock | rapid (record, reads); not stated | 12,824 | 0.260 → 0.079 (3.28) | 0.027 (0.021 to 0.033) | 0.026 (0.014 to 0.038) | +0.0008 |
| SRR30694000 | livestock | native barcoding (reads); ligation LSK112 (record); inconsistent record | 10,736 | 0.084 → 0.019 (4.44) | 0.007 (0.001 to 0.014) | 0.007 (-0.003 to 0.016) | -0.0001 |
| ERR16013022 | livestock | native barcoding (reads); not stated | 10,203 | 0.061 → 0.011 (5.53) | 0.006 (-0.002 to 0.022) | 0.007 (-0.020 to 0.034) | -0.0011 |
| SRR29786052 | wastewater | ligation (record, reads); R10.4.1 | 18,179 | 0.470 → 0.126 (3.74) | 0.076 (0.066 to 0.088) | 0.074 (0.049 to 0.099) | +0.0021 |

Reported next to the verdict, without a rule (amendment 06 Part C):

- **Both stages together:** 14 of 14 counting runs meet H1a, and 14 of 14 meet H2 (3 from stage 1, 11 from stage 2).
- **By class (stage 2):** human gut 6 of 6, livestock 4 of 4, wastewater 1 of 1.
- **By chemistry (stage 2):** R10.4.1 stated, 3 of 3; R9.4 stated, 1 of 1; chemistry not stated or record
  inconsistent, 7 of 7.

## What the result does and does not show

- **H1a is shown in every run.** Cutting the reads to an N50 of 3 kb lowers the naive share 2.0- to 7.7-fold.
- **H2 is informative in 5 of the 11 runs.** These have a net linked fraction of 0.076 to 0.153, and it changes by
  0.002 to 0.006 after cutting (1% to 6% of the value). In the other 6 runs the fraction is below 0.03, so "stable"
  there says only that a value near zero stays near zero.
- **A small downward shift.** In all 5 of these runs, and in the 3 counting runs of stage 1, the estimate after cutting
  is lower than the native one, by 0.5% to 6%. Each difference is far inside the tolerance and the intervals, but
  the direction is the same in 8 of 8 runs. The paper should state this as a small bias of the cut, not as exact
  invariance.
- **Wastewater contributes one counting run.** Five of the six accepted wastewater runs do not count, because sludge
  reads hold too few ARG (3 to 145 calls in 3 Gbp). The feasibility note of amendment 06 expected livestock to be the
  weak class; it was wastewater.

## Accepted runs that do not count

| Run | Project | N50 | ARG calls | Reason |
|---|---|---|---|---|
| ERR14656553 | PRJEB83983 | 13,505 | 145 | half-width at 3 kb is 0.12 (limit 0.10); 35 sides reach 2 kb after the cut |
| SRR27730983 | PRJNA1067729 | 12,804 | 41 | 4 sides reach 2 kb after the cut (30 needed) |
| ERR12160882 | PRJEB67571 | 13,910 | 65 | 18 sides reach 2 kb after the cut; half-width 0.12 |
| SRR11673973 | PRJNA629478 | 12,198 | 14 | 14 sides reach 2 kb in the native reads |
| DRR333437 | PRJDB12729 | 10,488 | 3 | 5 sides reach 2 kb in the native reads |

All five are wastewater. The four with a defined naive ratio have 1.85 to 2.38; the two with an estimate after the
cut are stable (ERR14656553 0.209 and 0.216; ERR12160882 0.144 and 0.123). They do not enter the verdict.

## The walk

25 candidates were tried and 16 accepted: human gut 6, livestock 4, wastewater 6. No fetch failed. No statistic was
computed on any run before the walk had ended and the kit rule had been applied
(`public/screen/p3x_status.tsv`, `p3x_manual_exclusions.tsv`, `p3x_project_kits.tsv`, `p3x_kit_probe.tsv`).

Candidates not accepted:

| Run | Class | Project | Subsample N50 | Reason |
|---|---|---|---|---|
| SRR36719965 | livestock | PRJNA1398240 (swine slurry) | 14,640 | kit rule: kit and chemistry both unknown |
| SRR36719968 | livestock | PRJNA1398240 | 12,719 | kit rule, as above |
| SRR36719970 | livestock | PRJNA1398240 | 13,562 | kit rule, as above |
| SRR38674067 | livestock | PRJNA1466632 (cattle rumen) | 11,753 | kit rule, as above |
| SRR38674069 | livestock | PRJNA1466632 | 11,751 | kit rule, as above |
| SRR24197343 | wastewater | PRJNA955434 (hospital wastewater) | 11,084 | kit rule, as above (ligation adapter at 1.7% of read starts, below the 5% threshold) |
| ERR14218439 | livestock | PRJEB83989 (cattle rumen) | 8,130 | N50 below 9 kb |
| ERR14218453 | livestock | PRJEB83989 | 7,714 | N50 below 9 kb |
| ERR14218448 | livestock | PRJEB83989 | 7,606 | N50 below 9 kb |

- **Kit rule.** For the six rejected runs the reads are trimmed or show no adapter or barcode at 5% of read starts,
  the archive records state neither kit nor chemistry, and no publication was found on PubMed or the web on
  2026-11-17. Each run was judged on its own reads.
- **Livestock stopped at 12 candidates tried** with 4 accepted, so stage 2 has 16 runs, not 18.

## Points a reader should know

1. **SRR29786052 (sewage influent).** Its record says the sample was enriched for bacteria with capture particles
   before DNA extraction. This is not PCR or whole-genome amplification and not an enrichment culture, and the project
   was in the frozen list, so the run stayed. Without it the count is 10 of 10 for both hypotheses.
2. **ERR12160882.** The submitter trimmed and filtered the reads (file name in the record). The run does not count.
3. **Chemistry is not stated for 6 of the 11 counting runs**, and one more record is inconsistent. The kit rule asks
   for kit or chemistry; the kit is known for all 11.
4. **Control strand in one run.** In SRR33225186, 2,462 of the 2,645 reads that Kraken2 assigns to *Escherichia coli*
   are 3.4–3.6 kb long, the length of the ONT lambda control strand (0.8% of the reads). They sit in the *E. coli*
   stratum and carry no call. The run has 9,527 ARG calls, so the effect on its estimate should be very small; it was
   not measured, because public reads are not stored. No other stage-2 run shows this peak.
5. **Runs were chosen for long reads after stage 1 was known.** This is why the stage-2 verdict is reported next to
   the stage-1 verdict and not in its place.

## Descriptive lines of the frozen script (not tested in stage 2)

- **H1b.** The script prints "not supported" (permutation p = 0.0535; within-class Spearman: gut 0.26, livestock 0.40,
  wastewater 0.71). Amendment 06 does not test H1b again, and the stage-1 result (supported, p = 0.0475) stands. Both
  values are close to 0.05, so the evidence for H1b is weak in either stage. Stage-2 runs span a narrow N50 range
  (10–19 kb), which limits what this line can show.
- **H3, pooled ARG:MGE within 2 kb.** Gut 0.124 (−0.097 to 0.345); livestock 0.023 (−0.034 to 0.081); wastewater
  0.079 (0.046 to 0.111); q 0.0045, 0.20 and 0.16. On the same read: gut 0.338, livestock 0.131, wastewater 0.295
  (q 0.0003, 0.035, 0.0003). Direction as in stage 1; the intervals are wide because each class has 4 to 6 runs.

## Deviations

None from the rules of amendment 06. Notes on timing and on the judgement calls above are in
`docs/preregistration_amendment_06_deviations.md`.
