# Gate G2 memo (amendment 03 v2)

2026-11-15. Frozen text: `docs/preregistration_amendment_03.md` (sha256 `125c206f…`); package commit `7b665ca`;
frozen scripts in `frozen/amendment03/`. Run 07:25–08:43 UTC. Verdict tables: `tables/G2/G2_verdict.tsv` and
`tables/G2/G2_g2b_settings.tsv`; rates with intervals: `tables/G2/G2_rates.tsv`.

## Verdict: G2 fails

| Part | Rule | Result | |
|---|---|---|---|
| G2a | ≥ 75% of the cases with E(N2-ev) ≥ 10 have \|ratio − 1\| ≤ 0.05, and the median ≤ 0.05 | 4 of 7 (57%); median 0.035 | fail |
| G2b | in each of 6 settings: two-sided rate within 3–6%, enrichment rate ≤ 6% | 3 of 6 settings pass | fail |

The pre-registered consequences (amendment 03, Part 3) apply and have been carried out:

- G2 does not stop the project.
- G2a fails, so **N2-ev is the primary read-level test** (exact under its own permutation null). N0-ev is reported as
  its smooth approximation.
- G2b fails in settings 0–2, so the rates and the cause are reported below. The paper limits the analytic test
  (N0-ev) to the settings that passed and uses N2-ev where it is exact.
- Read-level p values stay invalid at high coverage until P2 adds locus clustering.

## G2b: false-positive rates under no linkage

Two-sided rejection at α = 0.05 among the datasets with E ≥ 5 (3,920–4,000 per setting); 95% Wilson interval.

| Setting | Read lengths | A per genome, genomes | Mean E | Gated N0-ev | Mean O/E | Enrichment side | Pass |
|---|---|---|---|---|---|---|---|
| 0 | soil 279 | 24, 200 | 7.5 | 6.3% (5.6–7.1) | 0.80 | 1.0% | no |
| 1 | wastewater SRR32022521 | 24, 200 | 20.1 | 10.1% (9.2–11.1) | 0.86 | 0.8% | no |
| 2 | gut SRR23278673 | 24, 200 | 27.9 | 10.4% (9.4–11.3) | 0.90 | 1.0% | no |
| 3 | soil 279 | 6, 800 | 6.1 | 3.7% (3.2–4.3) | 0.97 | 3.7% | yes |
| 4 | wastewater | 6, 800 | 17.7 | 4.7% (4.0–5.4) | 1.00 | 4.3% | yes |
| 5 | gut | 6, 800 | 26.3 | 5.3% (4.7–6.0) | 1.01 | 4.9% | yes |

In settings 0–2, 92–97% of the false rejections are on the depletion side.

Reported, not gated (same datasets):

| Setting | N2-ev | N2-ev mean O/E | N0-ev without strata | Registered N0 |
|---|---|---|---|---|
| 0 | 5.9% (5.2–6.7) | 0.80 | 4.5% | 39.7% |
| 1 | 5.3% (4.6–6.0) | 0.93 | 5.6% | 30.8% |
| 2 | 4.8% (4.2–5.5) | 0.96 | 5.5% | 26.4% |
| 3 | 3.7% (3.2–4.4) | 0.97 | 4.5% | 40.6% |
| 4 | 4.4% (3.8–5.1) | 0.99 | 6.0% | 32.6% |
| 5 | 4.0% (3.4–4.6) | 1.00 | 7.9% | 24.9% |

- **The registered length-only null** rejects in 25–41% of datasets that contain no linkage, which confirms D30.
- **N0-ev without strata** is pushed toward enrichment where hosts differ (mean O/E 1.07–1.10 in settings 3–5),
  as expected without host strata.
- The soil settings did not bind at the 3% lower limit, as had been feared (3.7%).
- About 9.5% of loci were on ≥ 2 reads (0.1× depth), so pseudo-replication was small by design.

## Cause (diagnostic run after the verdict)

Scripts: `scripts/P1_15b_g2_diagnose.py`, `P1_15c_g2_summary.py`; tables `tables/G2/G2_diagnostic.tsv` and
`G2_presence_curve.tsv`. The diagnostic reran the first 1,000 G2b seeds of settings 0–2 and reproduced the gated
values exactly.

N0-ev keeps the positions of one class fixed and models the other class as genes placed independently in the gaps
it leaves. It fixes the class carried by more molecules.

- **Settings 0–2:** A sat on more molecules, so N0-ev fixed A (in 92–100% of datasets) and modelled B, which the
  generator places in islands.
- **Settings 3–5:** N0-ev fixed B and modelled A, which the generator places independently, so the model was right.

| Settings 0–2 | Two-sided rate | Mean O/E |
|---|---|---|
| Gated rule (A fixed, clustered B modelled) | 6.6%, 9.1%, 11.3% | 0.79, 0.86, 0.89 |
| Forced to fix B (independent A modelled) | 4.1%, 3.9%, 4.7% | 0.97, 1.01, 1.00 |
| B placed independently (gated rule) | 4.4%, 4.5%, 4.3% | 0.99, 1.01, 1.01 |

A clustered class violates the model in two ways:

1. **The presence curve has the wrong shape.** Calls that come in islands make the probability that a molecule
   carries B rise more slowly with room than the one-rate curve 1 − exp(−λx).
   - The fitted curve under-predicts at little room: observed / fitted 1.16–1.38 in the lowest deciles.
   - It over-predicts at a lot of room: 0.92–0.95 in the highest decile.
   - Molecules that carry A are longer, so their expected count is too high.
   - N2-ev bins on room instead of fitting a curve, so it removes this part.
2. **Molecules carrying A have B less often than other molecules with the same room** (observed / fitted
   0.63–0.91 at equal room). An island needs more free space than one gene, so the gene-level room overstates where
   B can sit.
   - Binning cannot remove this part, so N2-ev is also biased toward depletion here (mean O/E 0.80, 0.93, 0.96).
   - N2-ev stayed within the 3–6% band because the counts were small (E 7–28). With more molecules it would
     reject too often.

**G2a shows the same direction on real data.** In the three failing cases N0-ev expected 8–15% more co-occurrences
than N2-ev:

- SRR34068997 ARG:MGE, ratio 1.15;
- SRR23278673 ARG:MGE, ratio 1.08;
- SRR32022302 ARG:MGE, ratio 1.08.

The modelled class there is ARG. All three pairs are strongly enriched under both nulls (O/E 3.9–8.9), so no
conclusion changes.

## Real data are more clustered than the simulation

Clustering index (`scripts/P1_15d_clustering_index.py`, `tables/G2/G2_clustering_index.tsv`): the share of
carrier molecules with ≥ 2 calls of the class, divided by the share a one-rate Poisson model predicts. It is about 1
for genes placed independently. It also rises when density differs between hosts, so it does not by itself
separate physical clusters from host effects.

| Class | Simulated A | Simulated B islands | 8 G1b runs |
|---|---|---|---|
| ARG | 0.99 | | 7–165 |
| MRG | | | 6–134 where ≥ 900 molecules carry MRG |
| MGE | | 8.2 | 3–11 |

Every class in every run is clustered (index ≥ 3), and ARG and MRG mostly more so than the simulated islands
(8.2). In real data **both** classes break
the independence that N0-ev and N2-ev assume for the modelled class. So the rule "fix the clustered class" is not a
general fix.

What this means now:

- In ARG:MGE (MGE fixed, ARG modelled), read-level O/E is expected to lean toward depletion by an amount not yet
  known.
- Enrichment calls are conservative with respect to this bias, but p values are still too small because of
  pseudo-replication.
- Depletion calls ("ARG and MGE avoid each other") are not trustworthy at the read level until this is fixed.

## Consequence for `linkstat test` v0.1 (done; package commit `0cd23d4`)

- N2-ev is the primary null (`null_model = N2-ev`): exact p values, mid-p, TOST via the tilted cell laws, and the
  Mantel–Haenszel odds ratio with its interval. N0-ev (`n0ev_*`) and the registered null (`registered_*`) are
  reported alongside.
- New column `clustering_index_modelled` and the flag "modelled class clustered (O/E may lean toward depletion)"
  when the index exceeds 2. All real pairs tested so far carry it.

## Decision needed (yours: it changes a registered definition)

How to remove the clustering bias from the read-level test:

- **(a) Fold it into P2 (recommended).** Amendment 04 would treat runs of calls of one class within a set distance
  as one unit, with the unit's span as its excluded volume, for both the presence model and the room. A new
  simulation gate (G3) would use clustered A **and** B and host heterogeneity, with new seeds.
  - P2 already builds the locus machinery this needs, and read-level p values are not used for claims until P2
    anyway.
- **(b) Fix it now,** before P2, with amendment 04 and G3 on their own.
- **(c) Document it as a limitation only;** the paper relies on the net linked fraction and reports O/E as
  descriptive.

The net linked fraction (the effect size) is not affected through its anchors: dropping every anchor that touches
a B call moves it by ≤ 0.001 in all 16 real pairs (step 1.16 report). Whether clustering of A itself biases it
under no linkage has not been tested; G3 should include that check.
