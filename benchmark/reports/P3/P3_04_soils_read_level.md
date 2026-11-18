# P3 step 3.4b: read-level test on our 8 deep soils

2026-11-16. Our soils are not public data and are outside the registration; the same test and classes are used
(amendment 04; ARG and MRG at ≥ 80%, MGE = mobility core). Host strata: rule B (D56). Per-soil values only: no site
comparison is made here (vko8soil rule 8). Table: `tables/P3/P3_04_soil_reads.tsv` (rule B, rule A, and no strata).

## Result

**At read level our soils show no linkage beyond chance, and they have almost no power to show any.**

| Pair | Events within 2 kb, all 8 soils: observed / expected | Same molecule: observed / expected | Smallest p (2 kb; molecule) |
|---|---|---|---|
| ARG:MGE | 13 / 12.1 | 24 / 29.6 | 0.015; 0.051 |
| ARG:MRG | 0 / 0.12 | 0 / 0.16 | 0.97; 0.96 |
| MRG:MGE | 3 / 1.5 | 6 / 3.9 | 0.09; 0.12 |

- **ARG:MRG:** no read in any soil carries both within reach, against 0.1–0.2 expected by chance. Nothing can be
  estimated; the result is "no events", not "at chance".
- **ARG:MGE:** the net linked fraction within 2 kb is between −0.002 and +0.003 in the five soils where it can be
  estimated. One soil (307) has p = 0.015 (q = 0.045 across its three pairs) with an estimate of 0.003 and an interval
  of −0.03 to 0.04; one nominal result in eight soils is what chance gives.
- **Counts:** 83–1,357 ARG calls and 22–104 MRG calls per soil, on reads with an N50 of about 5 kb.
- **Host rule:** rules A and B and no strata give the same picture (estimates differ by at most 0.001). With rule B,
  19–67% of ARG-carrying reads are unresolved and there are 2–11 effective hosts.

## Reading

- This agrees with the published pilot (5 co-carrying reads among 2,418 ARG reads, at chance): the tool reproduces the
  qualitative read-level result with its own test.
- The read level cannot carry the soil story. The contig level (step 3.4a), where one contig spans what several short
  reads cannot, is the evidence to use for our soils, as the project plan expected.
- The Ridder against Pervomaisk comparison waits for the frozen soil catalogue and the reuse permission.
