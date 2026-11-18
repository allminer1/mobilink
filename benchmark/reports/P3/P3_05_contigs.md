# P3 step 3.4a: contig-level test on our soils

2026-11-17. Our soils are not public data and are outside the registration; the same classes and test are used
(amendment 04). Per-soil values and sums over all eight soils only: no site comparison is made (vko8soil rule 8). The
Ridder sum in part 1 reproduces a published value of one site.

## Result

**On contigs, antibiotic-resistance genes (ARG) and metal-resistance genes (MRG) share a contig about as often as
chance predicts. The published depletion (0.14-fold) comes from the null used, not from the data.**

- The published count is reproduced: 4 contigs carry both classes, against 28.6 expected by label permutation.
- On the same loci, linkstat expects 4.4 contigs and observes 4.
- In a simulation with no linkage, label permutation reports "depletion" in 67–100% of datasets once genes form
  clusters or contigs differ in gene density. linkstat with contig clumps rejects in 0–3%.
- On our eight soils with the registered classes, no contig carries an ARG and an MRG, and almost none is expected
  (0.15). With ARG and MRG at ≥ 60% identity, 96 ARG units share a contig with an MRG against 89 expected.
- Power is low everywhere: no soil can show equivalence to chance.

## 1. The published pilot's loci (Ridder, 4 samples)

Input: catalogue Rv1 of `ridder_study/` (801 ARG and 349 MRG loci on the metaFlye contigs of 278–281; read only).
Scripts `P3_05a_pilot_loci.py`, `P3_05b_pilot_test.py`; tables `tables/P3/P3_05_pilot.{nulls,linkstat}.tsv`.

| Chance model | Contigs with both classes: observed / expected | Ratio | p (depletion) |
|---|---|---|---|
| Label permutation on loci (the paper's primary null) | 4 / 28.6 (published 28.9) | 0.14 | < 0.0001 |
| Label permutation on units (same-class loci within 1 kb merged) | 4 / 14.8 | 0.27 | 0.0001 |
| Label permutation on units, within the MAG | 4 / 12.5 | 0.32 | 0.0001 |
| Length-weighted placement (the paper's second null) | 4 / 11.3 (published 11.2) | 0.35 | 0.009 |
| linkstat unit-level O/E (ARG positions fixed, MAG strata) | 4 / 4.4 | 0.91 (95% interval 0.25–2.34) | two-sided 1.0 |

- **linkstat same-contig test:** 5 of 787 ARG units have an MRG on their contig, against 7.2 expected (two-sided
  p = 0.35). Ratio of shares 0.81 (90% host-jackknife interval 0.41–1.60). By the H4 rule of amendment 05 this is
  "not determined": neither depletion nor equivalence is shown.
- **Within 10 kb:** 1 event against 1.4 expected.
- **Per sample (linkstat same-contig, observed / expected):** 278: 3 / 2.9; 279: 0 / 2.3; 280: 0 / 0.7; 281: 2 / 1.3.
  No sample has p ≤ 0.05. The label-permutation null gives p ≤ 0.0002 in 278, 279 and 280, as published.
- **Why the nulls differ:** 32 contigs carry two or more units; 18 carry only ARG, 10 only MRG and 4 both. Same-class
  genes share contigs beyond 1 kb (van clusters, cnr/ncc, mer with sil). A label permutation counts that as the
  classes avoiding each other. linkstat asks a different question: given where the ARG units are and how long the
  contigs are, is an MRG on the same contig less often than the MRG density predicts?
- **Host strata do little here:** only 9–29% of the loci lie on binned contigs.
- **Contig names in the pilot:** `p9_build_master.py` maps the MAG contigs back to the assembly names but keeps
  Bakta's renamed names for the unbinned sets. With the assembly names for both, the co-carrying contigs are the same
  4 and the resistance contigs are 1,075 instead of 1,074. The published count is unaffected. The pilot's placement
  nulls looked up unbinned contig lengths by the wrong names; with the right names the expectation is 11.3.

## 2. Both nulls where the truth is known

Script `P3_05e_null_sim.py`; table `tables/P3/P3_05_null_sim.summary.tsv`. Contigs and MAGs of sample 278; 270 A and
125 B genes placed with no dependence between the classes; 100 simulated assemblies per scenario.

| Scenario | Label permutation: observed / expected | Label permutation: "depletion" at 0.05 | linkstat same-contig, no clumps (A→B; B→A) | linkstat same-contig, contig clumps |
|---|---|---|---|---|
| Uniform positions | 0.99 | 7% | 5%; 2% | 3%; 0% |
| Gene clusters wider than 1 kb | 0.02 | 100% | 9%; 19% | 0%; 1% |
| Contigs differ in gene density | 0.40 | 76% | 9%; 5% | 1%; 3% |

- Label permutation is valid only when genes are placed uniformly. Merging loci into units does not repair it (100%
  and 67%).
- linkstat's same-contig test was too liberal on contigs when every unit counted as its own clump (up to 19% for a
  nominal 5%). Units on one contig share the event "this contig carries B".
- **Setting used from here on:** all A units of a contig are one clump (`contig_clumps`, passed as `loci`). The test
  is then conservative (0–3%). The estimates do not change. This goes into the package as an option in step 3.6.

## 3. Our eight soils, per assembly

Scripts `P3_05c_*` (annotation with the read pipeline on 70 kb windows) and `P3_05d_soil_contig_test.py`; tables
`tables/P3/P3_05_soil_contigs.{inputs,linkstat,nulls}.tsv`. Molecule = contig; stratum = Binette bin; contig clumps.
Two assemblies per soil (myloasm; Flye + medaka). Two class sets: registered, and ARG and MRG at ≥ 60% identity
(MGE unchanged), because the registered rule leaves few calls on soil contigs.

**Calls per assembly (myloasm, 0.55–2.19 Gb):** registered ARG 10–121, MRG 5–18, MGE 776–2,909; at ≥ 60% ARG 114–657,
MRG 327–1,069.

**Same contig, summed over the eight soils (observed / expected; soils with p ≤ 0.05):**

| Pair | Registered, myloasm | Registered, Flye | ≥ 60%, myloasm | ≥ 60%, Flye |
|---|---|---|---|---|
| ARG:MRG | 0 / 0.15; 0 | 0 / 0.01; 0 | 96 / 89.1; 0 | 50 / 51.4; 1 (above) |
| ARG:MGE | 24 / 20.3; 0 | 11 / 9.4; 0 | 117 / 123.2; 1 (below) | 62 / 65.2; 0 |
| MRG:MGE | 7 / 1.5; 1 (above) | 3 / 0.4; 1 (above) | 384 / 259.0; 6 (above) | 215 / 154.6; 3 (above) |

**ARG:MRG at ≥ 60%, per soil (myloasm; same contig):**

| Soil | Observed / expected | p | Ratio of shares (90% interval) | ARG contigs unbinned | Effective hosts |
|---|---|---|---|---|---|
| 278 | 32 / 33.9 | 0.60 | 1.09 (0.78–1.52) | 55% | 3.7 |
| 279 | 8 / 7.0 | 0.54 | 1.01 (0.53–1.93) | 49% | 5.3 |
| 280 | 10 / 13.8 | 0.11 | 0.99 (0.62–1.56) | 38% | 7.5 |
| 281 | 11 / 7.8 | 0.20 | 0.90 (0.20–3.99) | 70% | 2.1 |
| 304 | 10 / 10.0 | 0.99 | 1.13 (0.25–5.13) | 72% | 2.1 |
| 307 | 12 / 10.2 | 0.51 | 1.41 (0.42–4.79) | 74% | 2.1 |
| 308 | 7 / 3.2 | 0.05 | 2.55 (0.53–12.3) | 86% | 1.5 |
| 310 | 6 / 3.2 | 0.13 | 1.12 (0.26–4.93) | 65% | 2.4 |

- **ARG:MRG:** with the registered classes nothing can be estimated ("no events", with 0.15 expected). At ≥ 60% no
  soil shows depletion; one Flye assembly (308) is above chance (7 / 1.8, p = 0.008). No interval lies within
  0.8–1.25, so no soil is "at chance" by the H4 rule either: all are "not determined".
- **ARG:MGE:** at chance in both class sets. One soil is below chance at ≥ 60% (280: 11 / 16.9, p = 0.013).
- **MRG:MGE:** above chance at ≥ 60% in 6 of 8 myloasm and 3 of 8 Flye assemblies (ratio of shares 1.4–1.8 in most of them, 3.2 in soil 308). Within 2 kb
  the excess is 26 / 13.9 and no soil reaches p ≤ 0.05; within 10 kb it is 153 / 78.5. With the registered classes
  there are 7 events against 1.5. MRG at ≥ 60% includes many transport homologs, so this is a lead, not a result.
- **Label permutation on the same assemblies** reports "depletion" for every pair: ratios 0.13–0.41 on loci, in 8 of
  8 soils for ARG:MGE and MRG:MGE at ≥ 60%. That includes the pair that linkstat finds above chance.
- **Without host strata** the sums move little (ARG:MRG 96 / 76.2; MRG:MGE 384 / 220.5).

## Limits

- 38–86% of ARG-carrying contigs are unbinned and share one "unresolved" stratum; 1.5–7.6 effective hosts per soil,
  so the intervals are wide.
- "Same contig" depends on contig length. linkstat conditions on the contig lengths, but a fragmented assembly still
  hides links between genes more than a contig apart.
- The contig annotation here is the read pipeline (DIAMOND blastx on CARD, BacMet2 and mobileOG-db). The pilot used
  RGI, AMRFinderPlus and BacMet2 on predicted proteins, so its catalogue (part 1) is much larger than the registered
  classes (part 3).
- The Ridder against Pervomaisk comparison waits for the frozen soil catalogue and the reuse permission.
