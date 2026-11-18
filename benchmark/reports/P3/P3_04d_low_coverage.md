# P3 step 3.4d: read-level test on the 27 low-coverage soil libraries

2026-11-17. Libraries of our first soil study (PRJNA1436833). Our soils are not public data and are outside the
registration; the same test and classes are used (amendment 04; ARG and MRG at ≥ 80%, MGE = mobility core). Per-library
and pooled values only: no comparison between sites is made (vko8soil rule 8). Tables:
`tables/P3/P3_04d_s1.libraries.tsv` (one row per library) and `tables/P3/P3_04d_s1.linkstat.tsv` (per library and
pooled, two variants).

## Result

**The 27 libraries together (10.1 Gbp) cannot show linkage or its absence: no estimate is possible for any pair.**

| Pair | Events within 2 kb, pooled: observed / expected | Same read: observed / expected | p (2 kb; same read) |
|---|---|---|---|
| ARG:MGE | 0 / 0.36 | 1 / 1.25 | 0.73; 0.99 |
| ARG:MRG | 0 / 0.00 | 0 / 0.03 | 1.00; 0.97 |
| MRG:MGE | 0 / 0.14 | 0 / 0.18 | 0.87; 0.83 |

- **No estimate.** The net linked fraction is undefined for all three pairs, pooled and in each of the 27 libraries
  (81 library-by-pair tests). The method flags "no events" and "few effective events" in every case.
- **Why.** The pooled reads hold 360 ARG, 65 MRG and 16,519 MGE calls. Of 714 ARG sides, 213 reach 2 kb on the read.
  Chance alone predicts less than one ARG–MGE pair within 2 kb, so zero observed is not evidence of depletion.
- **The one co-occurrence.** One read (library EKR1_r1b11, soil 299) carries an ARG and an MGE, against 0.05 expected in
  that library (p = 0.049) and 1.25 expected in the pool. One nominal result in 81 tests is what chance gives.
- **Without host strata** the picture is the same (ARG:MGE 0 / 0.97 within 2 kb; 1 / 2.67 on the same read).
- **Agreement with the deep soils.** The 8 deep soils of step 3.4b also showed no linkage beyond chance with little
  power. The low-coverage libraries add no information beyond that.

For the paper this is a worked example of the power limit: at a read N50 of 5 kb, 10 Gbp of soil reads give 0.36
expected chance events for ARG:MGE, so soil libraries of this kind need contig-level analysis or far more data.

## Libraries

36 libraries were sequenced (run 1 and run 2, barcodes 01–18; `tables/P3/P3_04d_s1_samples.tsv`). 27 are analysed here.
Left out: the 8 soils that were sequenced again deep (analysed from the deep data in step 3.4b, so each soil counts
once) and the control library.

| Locality | Libraries | Gbp kept | ARG | MRG | MGE |
|---|---|---|---|---|---|
| Belousovka | 6 | 1.93 | 69 | 11 | 2,747 |
| Glubokoe | 2 | 1.23 | 54 | 9 | 1,897 |
| Pervomaisk | 2 | 0.95 | 35 | 6 | 1,527 |
| Ridder | 1 | 0.40 | 10 | 8 | 657 |
| Shemonaikha | 5 | 2.72 | 89 | 8 | 5,068 |
| TMK | 5 | 1.09 | 38 | 11 | 2,018 |
| Ust-Kamenogorsk | 6 | 1.80 | 65 | 12 | 2,605 |

Per library: 84–817 Mbp kept, read N50 4,620–5,833, 1–35 ARG calls and 0–8 MRG calls.

## Method

1. **Annotation.** The run script of the public runs, unchanged (`frozen/amendment05/P3_02_public_run.sh`), on all raw
   reads of each library (`scripts/P3_04d_s1_annotate.sh`). All 27 libraries finished without a failure.
2. **Control reads.** The libraries hold the ONT lambda control strand. A read is removed if it has an alignment of
   ≥ 200 bases to the lambda genome (NC_001416.1; `scripts/P3_04d_s1_lambda.sh`). This removes 50,744 reads (0.7–4.7%
   of the reads of a library; 186 Mbp in total), of which 1,708 are chimeras of a control strand and a soil molecule.
   18 gene calls (1 ARG, 17 MGE) sat on removed reads, all outside the control strand, on the soil part of a chimera.
3. **Test.** `scripts/P3_04d_s1_test.py`: loader and test of the frozen stage-1 script; per library, and pooled with the
   library as stratum, so anchors never leave their library and host. Host strata are Kraken2 species, as for the
   public soil runs. These libraries have no bins, so host rule B of the deep soils (D56) cannot be applied.
   55–64% of the reads of a library have no species label; 11% of ARG-carrying reads are in the unresolved stratum.

## Deviation from the design of D77

D77 planned to remove control reads by their Kraken2 taxon (genus Lambdavirus). That removed nothing: the Kraken2
database of the run script (standard, 16 GB) assigns the control reads to *Escherichia coli*. The check on the first
two libraries showed this before the test was run, and alignment to the lambda genome replaced the taxon rule (D78).
The test itself was run once, after the change.

Consequence for other steps: in the public runs and the deep soils, lambda control reads, if present, would sit in the
*E. coli* stratum. This was not checked here. Public runs with rapid or ligation kits rarely carry the control strand,
and the control strand got no call of the three gene classes in 50,744 reads here, so an effect on the registered results is unlikely but not excluded.
