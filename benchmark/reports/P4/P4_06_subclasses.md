# Phase 4, step 4.6: which ARG, and which MGE

Date: 2026-11-22. Descriptive only: no test is attached to these values (amendment 07, Part G; decision D112).
The sub-classes were fixed from the call counts before any linkage value was computed.

## What was done

- Sample: the 241 runs of the frozen analysis (86 projects, 7 environment classes).
- Measure: the ARG–MGE net linked fraction within 2 kb, as registered. Runs of a class are pooled; the interval is
  the delete-one-project jackknife. The functions are those of the frozen script (`frozen/amendment07/P4_03_confirm.py`).
- 16 sub-classes, each estimated in every class (112 cells):
  - 2 kinds of MGE, with all ARG: integration/excision genes and transfer genes (major mobileOG-db category).
  - 6 resistance mechanisms, with all MGE (CARD): efflux, inactivation, target alteration, target protection, target
    replacement, other.
  - 8 drug families, with all MGE (CARD drug classes merged): multidrug (two or more families), tetracycline, MLS
    (macrolide, lincosamide, streptogramin), aminoglycoside, beta-lactam, fluoroquinolone, peptide, other.
- Scripts: `scripts/P4_06a_subclass_counts.py`, `P4_06b_subclass_linkage.py`, `P4_06c_subclass_summary.py`.
- Tables: `tables/P4/sub/subclass.classes.tsv` (all cells), `subclass.mge_difference.tsv`,
  `subclass.composition.tsv`, `subclass.influence.tsv`; per cell with its jackknife replicates in
  `tables/P4/sub/linkage/<sub-class>/<class>.tsv`. Supplementary Table S11: `manuscript/supplement/Table_S11_subclasses.tsv`.

87 of 112 cells have an estimate. The other 25 have fewer than 30 sides that reach 2 kb, mostly in marine,
freshwater and soil. 10 of the 87 rest on fewer than 100 sides and are marked with * below.

## 1. Kind of MGE

| Class | Integration/excision | Transfer | Difference (95% interval) | Projects with a positive difference |
|---|---|---|---|---|
| Activated sludge | 0.098 | 0.001 | +0.097 (0.050 to 0.144) | 12 of 12 |
| Livestock | 0.099 | 0.031 | +0.068 (−0.011 to 0.148) | 8 of 8 |
| Marine | 0.060 | 0.001 | +0.058 (−0.002 to 0.118) | 10 of 10 |
| Wastewater | 0.044 | 0.003 | +0.042 (−0.033 to 0.116) | 9 of 9 |
| Human gut | 0.058 | 0.025 | +0.033 (−0.030 to 0.096) | 19 of 19 |
| Freshwater | 0.017 | −0.000 | +0.017 (−0.014 to 0.048) | 15 of 15 |
| Soil | 0.028 | 0.040 | −0.013 (−0.050 to 0.025) | 1 of 22 |

"Projects with a positive difference" counts the delete-one-project replicates, so it shows that no single project
sets the sign; it is not a count of independent confirmations.

- In six of seven classes ARG are linked to integration/excision genes (transposases, integrases, recombinases) and
  hardly to transfer genes (conjugation). The interval of the difference excludes zero only in activated sludge.
- Soil is the exception: the two values are similar, and both intervals include zero.
- The two values add up to the registered class estimate within 0.003 in every class (for example human gut
  0.058 + 0.025 against 0.085).

## 2. Resistance mechanism

Net linked fraction to any MGE within 2 kb.

| Mechanism | Human gut | Livestock | Wastewater | Activated sludge | Marine | Freshwater | Soil |
|---|---|---|---|---|---|---|---|
| Target replacement | 0.136 | 0.188 | 0.268 | 0.322 | – | – | 0.397 |
| Inactivation | 0.265 | 0.243 | 0.117 | 0.120 | 0.283* | 0.022 | 0.277 |
| Target alteration | 0.112 | 0.068 | 0.027 | 0.051 | – | – | −0.002 |
| Target protection | 0.046 | 0.069 | 0.106 | 0.033* | – | −0.002* | 0.003 |
| Efflux | 0.013 | 0.108 | 0.017 | 0.083 | 0.016 | 0.011 | 0.012 |
| Other | −0.006 | −0.007 | −0.012 | – | – | – | −0.004* |
| Share of the ARG units that are efflux | 42% | 12% | 49% | 25% | 73% | 62% | 33% |

- Target replacement genes (mostly *sul* and *dfr*) and inactivation genes are the linked ones: 0.12 to 0.40 where
  an estimate exists, except inactivation in freshwater (0.02).
- Efflux genes are at or near chance in five classes (0.011 to 0.017). They are the largest group of ARG calls in
  human gut, wastewater, marine and freshwater.
- So the class estimate is a mixture. A class with many efflux genes gets a low ARG–MGE value even when its
  acquired genes are linked as strongly as elsewhere. Marine shows this: 0.062 for all ARG, 0.28 for inactivation
  genes (49 sides), 0.016 for efflux genes, which are 73% of its ARG units.

## 3. Drug family

| Family | Human gut | Livestock | Wastewater | Activated sludge | Marine | Freshwater | Soil |
|---|---|---|---|---|---|---|---|
| Aminoglycoside | 0.213 | 0.278 | 0.089 | 0.074 | – | – | 0.305 |
| MLS | 0.195 | 0.148 | 0.172 | 0.057 | – | – | – |
| Beta-lactam | 0.181 | 0.190 | 0.089 | 0.171* | – | 0.016 | 0.019* |
| Other | 0.101 | 0.162 | 0.062 | 0.234 | – | 0.096* | 0.092 |
| Tetracycline | 0.059 | 0.065 | 0.057 | 0.022* | – | – | 0.007* |
| Fluoroquinolone | 0.025 | – | 0.011 | – | – | – | – |
| Multidrug | −0.001 | 0.070 | 0.002 | −0.001 | 0.006 | −0.003 | 0.002 |
| Peptide | −0.008 | −0.012 | −0.002 | – | – | – | −0.024* |

- Multidrug genes (mostly efflux pumps and their regulators) and peptide resistance genes are at chance or slightly
  below it in every class.
- Aminoglycoside, MLS and beta-lactam genes are linked at 0.06 to 0.28 in human gut, livestock, wastewater and
  activated sludge. Beta-lactam genes are near chance in freshwater and soil (0.02).

## 4. Single projects

Several of the high values depend on one project. Estimate, then the estimate without the most influential project:

| Sub-class | Class | Estimate | Without one project | Project | What it holds |
|---|---|---|---|---|---|
| Target replacement | Soil | 0.397 | 0.073 | PRJNA1104006 | 340 of 383 calls (*sul2*, *sul1*) |
| Aminoglycoside | Soil | 0.305 | 0.125 | PRJNA1171364 | – |
| Inactivation | Soil | 0.277 | 0.101 | PRJNA1171364 | 891 of 1,063 calls; 889 are one gene, APH(3')-Ia |
| Aminoglycoside | Human gut | 0.213 | 0.061 | PRJNA1377628 | – |
| Other families | Livestock | 0.162 | 0.055 | PRJNA1123871 | – |
| Aminoglycoside | Livestock | 0.278 | 0.183 | PRJNA1123871 | – |
| Efflux | Livestock | 0.108 | 0.042 | PRJNA1123871 | – |

- The three high soil values are each one project. One of them is one gene: APH(3')-Ia is a kanamycin resistance
  gene that is also a common marker on cloning vectors, so these reads may not be soil bacteria at all. This was not
  followed up; the registered analysis keeps the project.
- The wide jackknife intervals already show this (for example soil target replacement: −0.25 to 1.04).

## What this adds to the paper

1. The kind of MGE: where ARG are linked to MGE, the partner is an integration/excision gene, not a transfer gene.
   A 2 kb window sees transposons and integrons; it does not see whether the whole element sits on a conjugative
   plasmid.
2. The class value mixes intrinsic and acquired genes. Two environments can differ in the ARG–MGE value only because
   they differ in the share of efflux genes. Comparisons between environments should be made per mechanism or
   family, or with the composition stated.
3. With mechanism fixed, the classes look more alike than the overall estimates suggest (inactivation: 0.12 to 0.28
   in six classes), but the intervals are wide and several values rest on one project.
4. None of this is a registered result. It gives hypotheses for a later test, and it does not change the verdicts
   of the confirmatory analysis.

## Note on the run

The first complete run stopped at 82 of 112 cells: in marine, the only call of the mechanism "other" shares its ORF
with an MGE call, so no ARG unit was left and the estimator raised an error. The script now records this cell as
"no estimate"; the 82 finished cells were kept and the other 30 were computed. The check variant with all ARG and
all MGE equals the registered class estimate (marine, 0.0616). The package should give a clear message in this case
(listed for version 1.0).
