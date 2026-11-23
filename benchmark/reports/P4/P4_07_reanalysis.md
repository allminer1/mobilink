# Phase 4, steps 4.7 and 4.7b: reanalysis of a published comparison, and one sample with several kits

Status on 2026-11-22: complete (Parts 1 and 2, and step 4.7b).
All of this is exploratory (amendment 07, Part G). No registered test is attached. The designs were fixed before the
results were seen (decisions D114 and D115 in `docs/analysis_plan.md`).

## Source

Go and Shrestha 2026 (bioRxiv, doi:10.64898/2026.10.27.762393) annotated 40 public nanopore datasets from eight
environments. An ARG–MGE association is an ARG and an MGE call on the same read, with no distance limit and no chance
model. The preprint gives no ranking of environments. It reports that 1% of about 89 million reads carry an ARG and
that 4% of these also carry an MGE, and it shows that the samples group by environment.

Their read-level table is public (Zenodo, doi:10.5281/zenodo.23492757, `filtered_data.csv`, one row per ARG-carrying
read). It is in `external/go_shrestha_2026/` (md5 checked).

By the archive records:

- 14 of the 40 runs are amplified libraries: 9 by whole-genome amplification (wastewater treatment plants, swine
  wastewater) and 5 by random PCR (activated sludge).
- Each environment comes from one or two projects.
- None of the 40 runs is in our Phase 3 or Phase 4 samples.

## Part 1: their calls, and read length

Script `scripts/P4_07a_their_calls.py`; tables in `tables/P4/reanalysis/their_*.tsv`.

The table holds 3,443,876 ARG-carrying reads. 3.93% carry one of their MGE calls (insertion sequence, Tn3, Tn7 or
Tn554 marker, integron integrase). This matches the 4% of the preprint. 98.3% of the ARG reads have a plasmid
alignment, so the plasmid calls separate nothing and are kept apart.

**Within runs, the share rises steeply with read length.**

| Length of the ARG read | ARG reads | Share with an MGE call |
|---|---|---|
| under 1 kb | 56,108 | 0.0002 |
| 1–2 kb | 370,364 | 0.0044 |
| 2–4 kb | 893,492 | 0.0160 |
| 4–8 kb | 1,301,704 | 0.0387 |
| 8–16 kb | 758,500 | 0.0785 |
| over 16 kb | 63,708 | 0.1512 |

The log odds of an MGE call rise by 0.883 per doubling of read length (standard error 0.004; logistic regression
with one intercept per run, 40 runs). The rise is seen in every environment (`their_length_classes.tsv`).

**Across runs, the relation is weaker.** The Spearman correlation of the share with the mean read length of the run
is 0.36 (p = 0.02, 40 runs). With the median length of the ARG reads it is −0.01. Differences between projects are
larger than the length effect between runs.

**Cut to one length.** The ARG reads of every run were cut in silico to an N50 of 2 kb (Poisson breakpoints, 5
draws; a call stays on a piece when at least 80% of it is there).

| Environment | Runs | Projects | Amplified runs | Share, native | Share, cut to 2 kb | Rank, native | Rank, cut |
|---|---|---|---|---|---|---|---|
| Wastewater treatment plants | 5 | 1 | 5 | 0.201 | 0.074 | 1 | 1 |
| Hospital wastewater | 5 | 2 | 0 | 0.147 | 0.051 | 2 | 2 |
| Community sewage | 5 | 1 | 0 | 0.055 | 0.014 | 3 | 3 |
| Municipal wastewater | 6 | 1 | 0 | 0.029 | 0.006 | 4 | 6 |
| Swine wastewater | 4 | 1 | 4 | 0.027 | 0.007 | 5 | 4 |
| River | 5 | 2 | 0 | 0.023 | 0.003 | 6 | 8 |
| Activated sludge | 5 | 1 | 5 | 0.023 | 0.006 | 7 | 5 |
| Soil | 5 | 2 | 0 | 0.016 | 0.004 | 8 | 7 |

Shares are pooled over the reads of an environment.

- The share falls 2.7- to 7.5-fold in every environment.
- The three highest environments keep their order. The other five change places. The Spearman correlation between
  the two orders is 0.83.
- The pooled share of hospital wastewater rests on one run: the five runs range from 0.016 to 0.174, and the mean of
  the runs is 0.075.

**What Part 1 can and cannot say.**

- It can say that the published statistic depends on read length inside every run, and that its size is not
  comparable between datasets with different read lengths.
- It does not show that the order of the environments is an artefact: the top three keep their order after the cut.
- The cut uses only the ARG-carrying reads, because the table holds no other reads. A chance model needs all reads;
  that is Part 2.
- Environment, project and library type are confounded in this sample. Two of the three highest environments are one
  project each, and one of them is amplified throughout.

## Part 2: the 40 runs with our pipeline

The data pass ran from 02:37 to 10:34 UTC on 2026-11-22 (frozen run script of amendment 05: 3 Gbp of random blocks
per run, or the whole run; registered classes). All 40 runs were fetched; none failed. Scripts
`scripts/P4_07c_estimates.py` and `P4_07e_compare.py`; tables `tables/P4/reanalysis/ours.{groups,runs}.tsv` and
`compare.{env,stats}.tsv`. Runs of an environment are pooled.

| Environment | Runs | Projects | Read N50 (kb) | ARG units | Their share | Our naive share | Net linked fraction within 2 kb (95% interval) | Rank: theirs, naive, net |
|---|---|---|---|---|---|---|---|---|
| Wastewater treatment plants * | 5 | 1 | 7.8 | 15,378 | 0.201 | 0.358 | 0.051 (0.027 to 0.075) | 1, 3, 5 |
| Hospital wastewater | 5 | 2 | 3.1 | 5,315 | 0.147 | 0.419 | 0.280 (0.212 to 0.347) | 2, 1, 1 |
| Community sewage | 5 | 1 | 5.0 | 13,233 | 0.055 | 0.183 | 0.062 (−0.050 to 0.173) | 3, 5, 4 |
| Municipal wastewater | 6 | 1 | 4.9 | 136 | 0.029 | 0.104 | 0.020 (−0.011 to 0.051) | 4, 6, 6 |
| Swine wastewater * | 4 | 1 | 6.2 | 1,230 | 0.027 | 0.371 | 0.209 (0.160 to 0.257) | 5, 2, 3 |
| River | 5 | 2 | 6.4 | 176 | 0.023 | 0.036 | 0.007 (−0.064 to 0.078) | 6, 7, 7 |
| Activated sludge * | 5 | 1 | 5.3 | 138 | 0.023 | 0.221 | 0.221 (0.101 to 0.341) | 7, 4, 2 |
| Soil | 5 | 2 | 4.1 | 439 | 0.016 | 0.002 | 0.001 (−0.005 to 0.006) | 8, 8, 8 |

\* all runs are amplified libraries. "Their share": share of their ARG reads with one of their MGE calls (Part 1).
"Our naive share": share of our ARG-carrying reads with one of our MGE calls. The interval is the one of the package
(delete-one-host jackknife); the delete-one-run jackknife is in the table file.

**What the table shows.**

- **The published order changes with the annotation.** The order by their share and the order by our naive share
  agree only moderately (Spearman 0.62, p = 0.10, 8 environments; per run 0.56, p = 0.0001, 40 runs). Swine
  wastewater is fifth by their share and second by ours; activated sludge is seventh and fourth. The two annotations
  differ in the ARG rule (their alignment with at least 50% reference cover; ours at least 80% identity and 80%
  cover) and in what an MGE is (their insertion sequences, transposon markers and integrases; our mobileOG-db
  mobility core).
- **Chance correction changes the top of the order.** Wastewater treatment plants have the highest published share
  (0.20) and a naive share of 0.36 with our calls, but a net linked fraction of 0.051: most of the co-occurrence on
  the same read is what the density of MGE calls and the read length give by chance. They rank fifth by the net
  linked fraction. The naive share is 1.0 to 7.0 times the net linked fraction, depending on the environment.
- **With one annotation, the two measures give a similar order here.** Our naive share and the net linked fraction
  have a Spearman correlation of 0.88 (p = 0.004). The read N50 of the environments lies between 3.1 and 7.8 kb and
  is unrelated to either measure in this sample (Spearman 0.05 and −0.14). So in this dataset the read length is
  not what separates the environments; the annotation and the chance correction are.
- **The estimates are imprecise and rest on few data.** Four of the eight intervals exclude zero. Three
  environments have fewer than 200 ARG units in 3 Gbp per run (municipal wastewater, river, activated sludge;
  activated sludge has 54 sides that reach 2 kb). Only 18 of the 40 runs give an estimate on their own. Every
  environment is one or two projects, and three are amplified throughout, so no order of environments is supported
  by this sample under either measure.

**What Part 2 can and cannot say.**

- It supports two statements of the paper: the size of the naive share is not a measure of linkage (0.36 against
  0.05 in one environment), and a ranking from one or two projects per environment is not a ranking of environments.
- It does **not** show that the naive share and the net linked fraction order these eight environments differently.
  With our calls the two orders are close. The draft says so.
- The comparison with "their share" mixes two differences (annotation and subsample: their table covers whole runs,
  ours 3 Gbp per run). It was not designed to separate them.
- The partial result on 25 runs was seen on 2026-11-22 05:30 UTC while the script was tested; the design and the
  scripts were not changed afterwards.

## Step 4.7b: one sample, several kits

The data pass ran from 10:35 to 21:14 UTC on 2026-11-22 (81 runs fetched, 7 already held; none failed). Scripts
`scripts/P4_07c_estimates.py` (cells) and `P4_07d_kit_slopes.py`; tables `tables/P4/kits/kits.{groups,runs,slopes,
cells,profiles}.tsv`. A cell is one sample under one method, with replicate runs pooled. The slope is the change per
doubling of the read N50 of the cell, within samples (one intercept per sample, weighted by ARG reads or ARG units,
standard error clustered by sample).

| Set | Samples | Methods | Read N50 of the cells (kb) | Largest N50 ratio within a sample, median (range) |
|---|---|---|---|---|
| Stool, PRJNA1260441 | 10 (5 infants, 2 time points) | 2 extraction kits | 1.9–6.9 | 1.3 (1.0–1.9) |
| Rumen, PRJNA1117518 | 5 cows | 2 to 4 extraction or library kits | 1.2–7.5 | 2.0 (1.5–5.2) |
| Rumen, PRJNA1337107 | 4 cows | 4 extraction variants | 1.5–10.3 | 2.3 (1.6–6.8) |

**Slopes per doubling of read N50, within samples (90% interval).**

| Pair | Set | Cells / samples | Naive share | Net linked fraction within 2 kb | Cells / samples (net) |
|---|---|---|---|---|---|
| ARG–MGE | All sets | 48 / 19 | +0.016 (−0.025 to +0.058) | −0.015 (−0.048 to +0.017) | 41 / 16 |
| ARG–MGE | Stool | 20 / 10 | +0.009 (−0.081 to +0.100) | −0.032 (−0.094 to +0.030) | 20 / 10 |
| ARG–MGE | Rumen kits | 12 / 5 | +0.049 (+0.035 to +0.063) | +0.005 (2 samples only) | 5 / 2 |
| ARG–MGE | Rumen extraction | 16 / 4 | +0.002 (−0.065 to +0.069) | −0.003 (−0.024 to +0.017) | 16 / 4 |
| MRG–MGE | Stool | 20 / 10 | +0.042 (+0.022 to +0.062) | −0.002 (−0.013 to +0.009) | 20 / 10 |

MRG–MGE has too few MRG reads in the two rumen sets for a slope.

**What the table shows.**

- **ARG–MGE: no clear answer.** Over all sets the naive share rises by 0.016 per doubling, but the interval includes
  zero. It rises clearly in one set of three (rumen kits, +0.049). The net linked fraction has slopes near zero or
  negative, and none differs from zero; the interval over all sets (−0.048 to +0.017) is wider than the ±0.04 that
  the registered within-project analysis used as its bound.
- **MRG–MGE in stool: the expected pattern.** The naive share rises by 0.042 per doubling (p = 0.004) and the net
  linked fraction does not change (−0.002, interval −0.013 to +0.009).
- **The methods change more than the read length.** The genus profiles of one sample under two methods differ by a
  Bray–Curtis dissimilarity of 0.08 (stool), 0.15 (rumen extraction) and 0.35 (rumen kits) at the median; between
  samples under one method the values are 0.82, 0.24 and 0.45. In the rumen sets the method changes the recovered
  community almost as much as the animal does. The number of ARG reads also differs several-fold between methods
  (343 to 1,799 per method in the extraction set). So "same community, different read length" holds for the stool
  set only, and there the read N50 differs by a factor of 1.3 at the median.

**What step 4.7b can and cannot say.**

- It does not confirm, for ARG–MGE, that the naive share follows read length when the kit changes. One set shows it;
  the two others and the pooled estimate do not.
- It does not contradict it either: the intervals are wide, and the stool set has little spread in read length.
- For MRG–MGE in stool it shows the pattern that the in-silico cut shows.
- Public kit comparisons of this size are too small, and confounded with what the kit extracts, to replace the
  in-silico cut as evidence. The draft reports the result in these terms.

