# Phase 4, step 4.5: the registered cross-environment analysis

Date: 2026-11-21. Registration: `docs/preregistration_amendment_07.md` (frozen 2026-11-18). The analysis script
(`frozen/amendment07/P4_03_confirm.py`, package commit `86e3696`) was run once, from 14:06 to 21:19 UTC, by
`scripts/P4_03_run.sh`. Tables: `tables/P4/confirm/`; verdicts: `tables/P4/confirm/verdict.txt`
(sha256 9999aaa9b244…). Notes on the data pass: `docs/preregistration_amendment_07_deviations.md` (notes 1–7; no
deviation).

## Summary

- **H3a (primary) is not supported.** ARG are not linked to MGE more strongly in wastewater or in the human gut
  than in soil. Both differences are close to zero and far smaller than their standard errors.
- **Gate G5 is not met.** The registered fallback is the tier A paper with the Phase 4 results as they are. The
  decision is the user's.
- ARG sit next to MGE more often than chance in all 7 environment classes.
- The classes differ from each other (H3b supported), but no single pair of classes differs after correction.
- Within projects the naive share rises steeply with read length; the net linked fraction rises little, and its
  interval ends just outside the registered bound (H1c not determined for the net linked fraction).
- The differences between projects of one class are as large as the differences between classes. This is the main
  reason why the class contrasts are imprecise.
- No registered sensitivity analysis changes a verdict. Read-level filters change the class estimates little; the
  class definitions change their size.

## Sample

241 runs of 86 projects, all fixed from archive metadata before any read was fetched; 3 Gbp were subsampled from
each run. Of the 282 runs first selected, 41 were lost and no replacement run was accepted.

| Class | Runs | Projects | Effective projects¹ | Median read N50 (kb) |
|---|---|---|---|---|
| Human gut | 64 | 19 | 6.4 | 5.6 |
| Soil | 63 | 22 | 5.5 | 6.1 |
| Freshwater | 33 | 15 | 5.2 | 7.7 |
| Livestock | 26 | 8 | 3.4 | 6.7 |
| Activated sludge | 21 | 12 | 7.5 | 6.3 |
| Marine | 18 | 10 | 2.1 | 6.6 |
| Wastewater | 16 | 9 | 4.9 | 5.2 |

¹ By ARG units: the number of equally weighted projects that would give the same concentration. In every class a
few projects hold most of the ARG.

Rejected: 52 runs by the kit rule (19 amplified libraries; 33 with kit and chemistry unknown from reads, record and
publication) and 62 by fetch failure (58 of them are all the runs of one soil project that NCBI did not serve).

## Results

### Class estimates (ARG next to MGE, net linked fraction within 2 kb)

| Class | Class estimate | 95% interval (jackknife over projects) | Pooled test, q | Cut to 3 kb | Naive share |
|---|---|---|---|---|---|
| Livestock | 0.130 | 0.049 to 0.210 | 0.0003 | 0.123 | 0.331 |
| Activated sludge | 0.100 | 0.055 to 0.146 | 0.0003 | 0.101 | 0.222 |
| Human gut | 0.085 | −0.003 to 0.173 | 0.0003 | 0.072 | 0.221 |
| Soil | 0.069 | −0.020 to 0.158 | 0.0003 | 0.083 | 0.302 |
| Marine | 0.062 | −0.005 to 0.129 | 0.0002 | 0.070 | 0.109 |
| Wastewater | 0.049 | −0.047 to 0.145 | 0.0003 | 0.066 | 0.225 |
| Freshwater | 0.016 | −0.014 to 0.046 | 0.019 | 0.011 | 0.038 |

- The pooled test asks whether the ARG of these runs are linked to MGE above chance. It is significant in all 7
  classes.
- The jackknife interval asks a different question: how well the runs describe the class when projects are the
  sampling unit. It is wide in five classes and includes zero there.
- Cutting all reads to 3 kb changes the class estimates by 0.000 to 0.017.
- The order of the classes by the naive share differs from the order by the net linked fraction. Soil is second by
  the naive share and fourth by the net linked fraction; wastewater is third and sixth.

### H3a (primary): not supported

| Contrast | Difference | SE | t (df) | Holm p (one-sided) |
|---|---|---|---|---|
| Wastewater − soil | −0.020 | 0.060 | −0.33 (23.8) | 0.79 |
| Human gut − soil | +0.016 | 0.060 | 0.26 (38.9) | 0.79 |

The intervals are wide (about ±0.12), so the data do not show that the classes are equal either. Three projects
drive the estimates:

| Project left out | Class | Share of the class's ARG units | Class estimate without it |
|---|---|---|---|
| PRJNA1104006 | Soil | 31% | 0.030 (from 0.069) |
| PRJEB76456 | Wastewater | 35% | 0.093 (from 0.049) |
| PRJNA1112760 | Human gut | 33% | 0.121 (from 0.085) |

Leaving one project out changes the sign of both contrasts: wastewater − soil runs from −0.034 to +0.024, human
gut − soil from −0.007 to +0.054.

### H3b (secondary): supported

The classes are not all equal (Q = 18.1, 6 df, p = 0.006). The largest differences are activated sludge −
freshwater (0.084; Holm p = 0.065) and livestock − freshwater (0.114; Holm p = 0.25). No pairwise difference is
significant after the Holm correction over the 21 pairs.

### H1c (secondary): not determined for the net linked fraction

Within projects (human gut, livestock, wastewater, activated sludge; runs with at least 50 ARG-carrying reads):

| Outcome | Slope per doubling of read N50 | Interval | Runs / projects | Rule |
|---|---|---|---|---|
| Naive share | +0.117 | SE 0.019; one-sided p = 1.2 × 10⁻⁶ | 96 / 25 | positive: met |
| Net linked fraction | +0.016 | 90%: −0.011 to +0.044 | 91 / 24 | within ±0.04: not met (upper end 0.044) |

The naive share rises seven times faster than the net linked fraction, but the registered bound for the net
linked fraction is missed by 0.004.

Descriptive, from the per-run cut columns: 25 runs met the counting conditions of amendment 05 (enough ARG and a
long enough read N50). In all 25 the naive share fell at least 1.5-fold after the cut to 3 kb while the net linked
fraction stayed within its tolerance (median absolute change 0.002).

### H4b (secondary): ARG and MRG on the same read

| Class | Observed / expected reads | Share ratio (90%) | q | Outcome |
|---|---|---|---|---|
| Activated sludge | 10 / 3.3 | 3.01 (1.98–4.59) | 0.043 | above chance |
| Human gut | 2,233 / 1,590 | 1.48 (1.22–1.79) | 0.18 | not determined |
| Wastewater | 1,227 / 1,054 | 1.48 (1.35–1.63) | 0.30 | not determined |
| Livestock | 80 / 62.8 | 1.49 (1.12–1.99) | 0.27 | not determined |
| Soil | 14 / 24.7 | 0.98 (0.46–2.08) | 0.20 | not determined |
| Marine | 3 / 0.7 | 7.22 (1.39–37.7) | 0.18 | not determined |
| Freshwater | 2 / 0.5 | 2.89 (0.49–17.2) | 0.20 | not determined |

### Other descriptive results

- **Per run.** 136 of 241 runs had enough ARG for an estimate; 79 were above chance at q ≤ 0.05 and 6 estimates
  were below zero. Soil, freshwater and marine runs often hold too few ARG in 3 Gbp (median 30, 12 and 21 calls).
- **Per project.** In the human gut the project estimates run from 0.004 to 0.264 (19 projects, all positive); in
  wastewater from 0.012 to 0.158 (9 projects). In soil only 8 of 22 projects have an estimate (median 0.007; one
  project at 0.373).
- **Equal project weights** (random-effects mean of the project estimates; a registered sensitivity analysis only
  where every project has an estimate): wastewater 0.100 (0.058–0.142), human gut 0.079 (0.053–0.104), livestock
  0.076 (0.021–0.131). The pooled and the equally weighted estimates differ most in wastewater, where one
  long-read project (PRJEB76456, N50 14 kb, estimate 0.012) holds a third of the ARG.
- **MRG next to MGE**: above chance by the same-read test in all 7 classes (q ≤ 0.013).

## Sensitivity analyses that change the sample (Part E; `tables/P4/sens/sens_samples.*.tsv`)

Run on 2026-11-21 after the verdicts (`scripts/P4_04_sens_samples.py`). None changes a verdict.

| Scenario | Runs left out | Soil estimate | Wastewater − soil | Human gut − soil |
|---|---|---|---|---|
| Registered | none | 0.069 | −0.020 (SE 0.060) | +0.016 (SE 0.060) |
| Runs public before 2021 left out | 12 (activated sludge 7, human gut 3, marine 1, soil 1) | 0.068 | −0.019 (SE 0.060) | +0.018 (SE 0.061) |
| Danish consortium left out¹ | 10 soil runs of 3 projects | 0.104 | −0.054 (SE 0.056) | −0.019 (SE 0.056) |
| Runs with megabase reads left out² | 4 soil runs | 0.075 | −0.026 (SE 0.061) | +0.010 (SE 0.061) |

- Without the runs before 2021, activated sludge stays at 0.101, human gut at 0.086 and marine at 0.061.
- The Danish soils have little ARG–MGE linkage; without them the soil estimate is as high as that of activated
  sludge, and both contrasts are negative.

¹ The amendment does not name the projects. They were read as the four projects with site ids of the Danish
national survey (PRJEB104930, PRJEB114667, PRJEB85558, PRJEB125553; the last has no accepted run). This reading was
fixed after the confirmatory results were known.
² Exploratory (note 4 of the deviations file): runs in which reads over 200 kb hold more than 1% of the bases.

## Sensitivity analyses that change the reads or the classes (Part E; `tables/P4/sens/sens_variants.*.tsv`)

Run on 2026-11-21 and 22 after the verdicts (`scripts/P4_04b_sens_variants.py`; class estimate and delete-one-project
jackknife by the functions of the frozen script). The script reproduces the registered estimates exactly when
nothing is changed. H3a is not supported in any variant.

| Variant | Soil | Wastewater | Human gut | Wastewater − soil | Human gut − soil | Holm p | A single project changes a sign |
|---|---|---|---|---|---|---|---|
| Registered | 0.069 | 0.049 | 0.085 | −0.020 (SE 0.060) | +0.016 (SE 0.060) | 0.79, 0.79 | Yes |
| Chimeric reads removed | 0.064 | 0.028 | 0.083 | −0.036 (SE 0.053) | +0.019 (SE 0.062) | 0.76, 0.76 | Yes |
| Control reads removed | 0.069 | 0.049 | 0.085 | −0.020 (SE 0.060) | +0.016 (SE 0.060) | 0.79, 0.79 | Yes |
| Genus instead of species strata | 0.081 | 0.052 | 0.087 | −0.029 (SE 0.064) | +0.006 (SE 0.063) | 0.93, 0.93 | Yes |
| MGE at 80% identity | 0.065 | 0.045 | 0.081 | −0.020 (SE 0.052) | +0.015 (SE 0.054) | 0.78, 0.78 | Yes |
| ARG and MRG at 60% identity | 0.014 | 0.033 | 0.057 | +0.019 (SE 0.022) | +0.043 (SE 0.029) | 0.20, 0.16 | No |
| MGE as the broad classes of D3 | 0.077 | 0.127 | 0.219 | +0.050 (SE 0.105) | +0.142 (SE 0.077) | 0.32, 0.074 | No |

- **Read-level filters change little.** Removing chimeric reads, removing control reads, genus-level strata and
  the stricter MGE threshold move every class estimate by 0.03 or less. The largest changes are activated sludge
  with MGE at 80% (0.100 to 0.071) and wastewater without chimeric reads (0.049 to 0.028).
- **The class definitions change the size of the estimates.** With ARG and MRG at 60% identity there are 1.5 to 10
  times more ARG, and all class estimates fall (soil 0.069 to 0.014; activated sludge 0.100 to 0.031). With the
  broad MGE classes the estimates rise in the host-associated and engineered classes (livestock 0.273, human gut
  0.219, activated sludge 0.151, wastewater 0.127) and change little in soil, marine and freshwater (by 0.015 or less).
- **Under both alternative definitions the two contrasts are positive and no single project changes a sign**, but
  neither is significant after the Holm correction. These variants are sensitivity analyses; the registered test
  used the registered classes.
- Chimera filter: reads flagged from the all-against-all overlaps or with an adapter inside the read; a median of
  1.0% of the reads of a run was removed (largest 12.8%). Only reads with overlaps can be judged (median 6.7% of
  the reads of a run); 9.4% of the judged reads were flagged (median over runs).

### Control reads (registered descriptive item; `tables/P4/sens/control_screen.tsv`)

37 of 241 runs (23 projects, all 7 classes) carry suspected control-strand reads: reads labelled *E. coli* by
Kraken2 with a length of 3,400 to 3,619 bases, in runs flagged by the rule of D92 (34 runs) or by an excess of
this length window among the reads labelled *E. coli* (37 runs). 192,350 of 186.0 million reads were removed (at
most 2.8% of a run). No class estimate changed by more than 0.0005.

## Gate G5

| Condition (approved 2026-11-18) | Result |
|---|---|
| H3a is supported | **No** |
| At least 5 classes are analysed | Yes (7) |
| No single project changes the sign of either H3a contrast | **No** (both signs change) |

G5 is not met. The registered consequence: "the fallback is the tier A paper with the Phase 4 results as they
are."

**Decision of the user, 2026-11-22 (D113):** the content of the fallback is kept (tool paper with the 241-run survey
and the results as they are), and it is submitted to a tier B journal first (Microbiome as a Methodology article,
then iMeta), then to the tier A journals. The manuscript reports H3a as not supported and states that G5 was not met.

## What the results can carry in a paper

1. The tool finds ARG–MGE linkage above chance in every environment class, in 241 public runs of 86 projects.
2. The naive share ranks the environments differently from the net linked fraction and depends on read length
   within projects; the net linked fraction depends on it much less.
3. A ranking of environments by ARG mobility is not supported by these data. Projects of one environment differ as
   much as environments do, and a few projects hold most of the ARG. This is a limit for any cross-study
   comparison of this kind, including those made with the naive share.
4. The size of the net linked fraction depends on how the classes are defined (identity threshold for ARG; narrow
   or broad MGE classes). Values are comparable between studies only under one definition.

## Still to do in Phase 4

- Part E is complete.
- Steps 4.6 (running), 4.7, 4.7b and 4.9; step 4.8 is dropped (D113).
