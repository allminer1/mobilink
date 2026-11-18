# Amendment 01 to the pre-registration of H1 to H4

Date: 2026-11-14. It adds to `docs/preregistration_H1-H4.md` (sha256 `b1212d828c6cf387825adeec48de7263aa0650a82ac429a33291f33f5cbfc8b6`) and does not change it. The primary analysis, thresholds, window, tests and fixed sensitivity analyses stay exactly as registered.

## What is added

**One exploratory sensitivity analysis for the MGE class.** The registered MGE class has five mobileOG gene families removed: recA, recX, copR, clpP and nth.
- It is applied at the primary thresholds (ARG 80/80, MRG 80/80, MGE 60/80).
- It covers the pairs ARG:MGE (all ARGs and without the intrinsic templates) and MRG:MGE.
- It is reported next to the registered result, never in place of it.
- Script: `scripts/P0_05b_mge_sensitivity.py`.

## Why, and when it was decided

It was decided **after** the P0 own-soil results had been seen (`reports/P0/step0.5_ridder_report.md`, section 4). In the 8 own soils, most non-intrinsic ARG–MGE and MRG–MGE co-occurrences paired a chromosomal gene with a mobileOG homolog of a housekeeping protein:
- rsmA with recA, at about 3.6 kb;
- pstB and vanX with copR-family response regulators.

These neighbours are conserved on chromosomes and say nothing about mobility. The five families were chosen from those results, so the analysis is exploratory and labelled as such.

## What is not decided here

A data-independent "mobility core" MGE definition is left to the start of P1, before any confirmatory test. It is to be built from rules on gene family and evidence level, not from our hits. If adopted, it will be registered as amendment 02.
