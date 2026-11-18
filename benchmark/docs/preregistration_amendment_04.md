# Amendment 04: the net linked test, units, and gate G3

Status: **approved by the user on 2026-11-16 ("approve") and frozen the same day, before any G3 statistic.** No G3 seed (90,000,000 and above) had been used.
- **Approved decisions** (marked **[approved 2026-11-16]** below):
  - the margin rule for net_w (2.3);
  - how H2–H4 read with the new test (2.6);
  - the density clip at ±2.5σ (G3a generator);
  - the pass rules R1–R4 (G3a).
- **Still open:** the read-to-bin rule for our soils (2.1). G3 does not use it. It is decided before P3.
- **History:** draft v1 was approved in principle (D36) before the feasibility check. Draft v2 (`docs/preregistration_amendment_04_draft_v2.md`) is the text the user approved. This file is that text with the decisions recorded and the freeze details added (status, run time on 32 processes, the G3b mix seed). Nothing else changed.
- **Frozen code:**
  - the package at commit `100c199fabe3c7b9a2de1867aad07f093aacc3da` (`frozen/amendment04/PACKAGE_COMMIT`);
  - the G3 scripts in `frozen/amendment04/`, the three read-length profiles and the driver, with sha256 in `frozen/amendment04/SHA256SUMS`;
  - `scripts/G3a_run.sh` checks this file's sha256, the frozen checksums, the package commit and that the package source is clean. It then runs G3a and the G3a and G3b summaries.
- The sha256 of this file is in `docs/preregistration_amendment_04.sha256`.
- The registered pre-registration (`docs/preregistration_H1-H4.md`) and amendments 01–03 are unchanged as texts. This amendment replaces the registered read-level test (Part 2) and adds gate G3.

It follows from:
- the G2 failure (`reports/P1/G2_memo.md`);
- D35: fold the clustering fix into P2;
- D39: the registered read-level test becomes the net linked test;
- D40: host labels from Kraken2 for public runs and MAGs for our soils.

The feasibility check (plan step 2.7) was done on development seeds only. The other values were set by the feasibility runs, under the delegation, and are logged in `docs/analysis_plan.md` (D41–D45). The freeze is D46.

## 0. Why the read-level test changes

- **O/E against an exposure model is fragile at high counts.** Every small error in the model of where a gene could sit on a read becomes a false positive once E is in the hundreds. At 5× depth, unit-room model error gave O/E 1.11 (soil read lengths) and 10–25% false positives (P2 step 2.4).
- **The net linked fraction does not need that model.** It compares the real A units with random anchors placed on the same host's molecules, so read lengths and excluded volume enter empirically.
- **Two kinds of dependence must be handled:**
  - **Hosts.** A-rich genomes are also B-rich, which gives an association without linkage. Anchors placed within host strata remove it.
  - **Reads that cover one locus.** They see the same neighbouring genes. Without accounting for this, the test rejects 41% under no linkage.
- **Inference that was tried and failed** (D42, D44):
  - a t reference on the effective number of hosts;
  - a sign flip on jackknife pseudo-values;
  - score tests with observed variance, by host or by locus cluster.

  All of them rejected 6–15% under no linkage, because a few hosts or a few rare events carry the variance.

## Part 1. Units (unchanged from draft v1, except where noted)

- **Unit.** On a molecule, A calls whose gaps are ≤ d = 1,000 bp and that have no B call between them form one unit. It spans from the first start to the last end. Sensitivity: d = 500 and 2,000 bp.
- The net linked test measures A as units: under the null, an island of A genes is placed at random as a whole.
- **Room formulas for O/E (draft v1 Part 1)** stay as written. N2-ev-u and N0-ev-u are now **descriptive** (reported with their flags); they are no longer the registered test.

## Part 2. The registered read-level test: the net linked test

### 2.1 Host strata

- **Public runs:** Kraken2 2.17.2 with Standard-16 (2024-12-28) on all reads.
  - `--confidence 0`, with at least 10 k-mers that hit the database per read.
  - The stratum is the species. A read classified only at genus or family keeps that taxon as a coarser stratum.
  - Reads that are unclassified, classified above family, or below 10 hit k-mers form one stratum, "unresolved", which is flagged.
  - Why confidence 0: Kraken2's confidence divides by all k-mers, and at about 90% read accuracy no read reaches 0.1. In the gut run, 1.5% of reads were classified at 0.1 and 69% at 0 (95% of carriers).
  - Sensitivity: genus strata; confidence 0.1 where at least half of the carriers are classified at 0.1.
- **Our soils:** reads assigned to their MAG or bin.
  - **[open; decided before P3, not used by G3]** the read-to-bin rule. The default is the bin of the contig that the read maps to (minimap2, mapq ≥ 20), else unresolved.

### 2.2 Anchors

- Each host stratum gets 40 anchors per A unit in it, with at least 20,000 anchors in all.
- They are placed on the stratum's molecules with weight (L − length)+, with lengths drawn from the A unit spans.
- An anchor that touches any B call is redrawn.

### 2.3 Effect sizes (descriptive, with 95% intervals)

- **net_w = 1 − (1 − F_A(w)) / (1 − F_R(w)), at w = 2 kb.** F is the weighted Kaplan–Meier CDF of the far-end distance to the nearest B (the frozen Part B estimator, on units).
- **net_mol = 1 − (1 − M_A) / (1 − M_R), same molecule.** M is the weighted share of A units (or anchors) on a molecule that carries a B call.
- **Margin rule, for both** **[approved 2026-11-16; new for net_w, where it changes how the frozen Part B estimator is applied to units]**: only units ending ≥ d from both molecule ends count, and anchors get the same margin. The weight is 1 / p_whole(span + 2d).
  - Why: a unit cut by a molecule end looks complete but sits on a read chosen by its span, so a partly covered island is counted as if it were a whole unit. With A in islands this biased net_mol by −0.01, and net_w under linkage by about −0.006 in soil (0.026 against a truth of 0.032; run 6). With the rule both match the truth (probe v7).
- **Intervals:** delete-one-host jackknife SE, estimate ± t(G_eff − 1) × SE, where G_eff is the effective number of hosts by A sides. Without strata: a bootstrap over molecules.

### 2.4 Unit loci (which reads see the same locus)

- minimap2 `-x ava-ont` overlaps among the carrier reads, with block ≥ 1,000 bp and nmatch / block ≥ 0.25. The density rule stops homologous genes at different loci from joining: median density 0.92 within a locus and 0.15 between loci (soil 279).
- An A unit that lies ≥ 80% inside an overlap is mapped to the other read. It joins the unit there if the two overlap by ≥ 50% of the shorter.
- A join through a reverse-strand overlap swaps left and right.
- A locus whose joins disagree on orientation has its two sides merged into one clump (net_w), or one shared draw for all its occurrences (same molecule). This can only overstate the variance.
- In simulation, units join by their genome intervals (the same rule, with exact overlaps).

### 2.5 Tests (the registered p values)

- **Sides.** For each A-unit side, the room is the distance to the molecule end on that side. The outcome is y = 1 if the far end of the nearest B on that side lies within τ = min(w, room).
- **Chance.** F_h is the weighted Kaplan–Meier CDF of the anchors in the side's host stratum.
- **Clumps.** The sides that see the same genome side of one unit locus share one distance D ~ F_h, so:
  - E_k = Σ F_h(τ_i);
  - Var_k = Σ_ij F(min(τ_i, τ_j)) − E_k².
- **net_w test.** z = (O − E) / √(Σ_k Var_k), with clump = unit locus × oriented side. It uses all A-unit sides, not only the margin-rule units, which keeps its power. Selecting sides by their position on the read does not depend on B, so the test stays valid.
- **Same-molecule test.** One outcome per unit occurrence: a B on either side within the room. Under the null, P(no B) = (1 − F(room_L)) (1 − F(room_R)). The clump is the unit locus. Left and right are genome sides (read sides swapped for reverse-strand joins), and the covariance follows from the shared D_L and D_R.
- **Reference: the exact null distribution of O under this clump model**, by Monte Carlo with 9,999 draws. Each clump draws its shared distances from F_h. Mid-p values: p_up = (#{O* > O} + (#{O* = O} + 1) / 2) / 10,000, p_low likewise; two-sided p = min(1, 2 min(p_up, p_low)); the enrichment p is p_up. The smallest two-sided p is 1 × 10⁻⁴.
  - Why not the normal reference: at few events O is skewed and discrete. In run 6 the normal net_w test rejected 4.7% on the enrichment side (at 2.5%) when E²/V < 5, and almost never on the depletion side.
  - z = (O − E) / √V and the effective events E² / V are reported. Below 10 events the result is flagged, not excluded.
- **Use.** The same-molecule test answers the co-occurrence question (it replaces the O/E test in H2–H4). The net_w test answers local linkage. Both are reported for every pair, with Benjamini–Hochberg across pairs within each level and each test (as registered).
- **Flags** (reported with every result):
  - dense classes (A + B calls cover > 5% of all bases);
  - no host strata;
  - no locus clumps;
  - locus conflicts (> 10% of units);
  - fewer than 30 hosts with A units;
  - fewer than 30 sides reach w;
  - share of carriers in the unresolved stratum.

### 2.6 How H2–H4 read with the new test **[approved 2026-11-16]**

The registration states H3 and H4 in O/E terms, with TOST on log(O/E) and bounds ±log 1.25. The reading:
- **"O/E" means the anchor-based ratio M_A / M_R on the same molecule** (margin rule). It is an O/E with the anchors as the chance model: 1 under no linkage, at any read length.
- **H3 (enriched):** the same-molecule test, enrichment side, and net_w > 0. The between-environment comparison is planned in P4.
- **H4 (at or below chance, two-sided):** the same-molecule test two-sided. "At chance" is TOST on log(M_A / M_R) with the registered bounds ±log 1.25, on the 90% host-jackknife interval.
- **H2:** unchanged in wording. The test of "O/E stays at 1" becomes the same-molecule test.

## Part 3. Gate G3

### G3a: calibration in simulation (G3 seeds, 90,000,000+)

- **Generator** (`linkstat simulate`):
  - genome strata;
  - B: 5 islands + 12 singletons per 4 Mb genome;
  - A independent (24) or in islands (12 + 3 islands);
  - host heterogeneity: per-genome density log-normal with σ = 1, **clipped at ±2.5σ** **[approved 2026-11-16]**;
  - abundance log-normal with σ = 1; chimeras off.
- **Why the density clip.** Without it, 2,000 genomes reach about 30× the base gene density. In genomes that crowded, the simulator's non-overlapping placement packs A genes against B genes, a hard-core contact effect. It creates association under no linkage: community truth up to +0.23 in one dataset, and 7.6% rejection in the hardest settings. With the clip it is 5.2%.
  - Real genomes carry ARGs and MGEs among thousands of other genes, so this packing has no direct real counterpart. The "dense classes" flag covers real data.
- **Settings:** 3 read-length profiles (soil 279, wastewater SRR32022521, gut SRR23278673) × {A independent, A in islands} × depth, giving 18 settings:
  - 0.1× with 2,000 hosts;
  - 5× with 50 hosts;
  - 20× with 50 hosts.
- **Datasets per setting:** 500 under no linkage, and 200 each at linked shares 0.03 and 0.1. Run time is about 5 h on 16 processes, or about 2.5 h on the 32 used (the number of processes does not change any result).
- **Seeds:** 90,000,000 + 1,000 × setting + {0, 500, 700} for linked share {0, 0.03, 0.1} + dataset index. Each dataset has its own seed (90,000,000–90,017,899). All development runs used seeds below 90,000,000.
- **Pass rules, for each of the two tests** **[approved 2026-11-16; they replace the draft rule "3–6% in every setting"]**:
  - R1 size: pooled two-sided rejection within 2.5–6%, and ≤ 8% in every setting;
  - R2 enrichment side (one-sided 2.5% level): pooled ≤ 3.0%, and ≤ 5.0% in every setting;
  - R3 bias: |mean net_w| ≤ 0.01 and |mean net_mol| ≤ 0.01 under no linkage, every setting;
  - R4 coverage: the net_w 95% interval covers the community's truth in ≥ 90% of datasets, every setting, under no linkage and at each linked share.
  - Why not "3–6% in every setting": with 500 datasets per setting, an exactly calibrated test lands within 3–6% in all 18 settings only 6% of the time.
  - Why the R1 floor is 2.5%: the net_w test is conservative (3.3–3.7% in runs 6 and 7; below). That is valid and costs power. The floor still catches a broken test, for example one whose variance is overstated twofold (size below 1%).
  - What the new rules do, for each test: an exactly calibrated test fails R1 in 3% of gates and R2 in 1%. A setting whose true size is 10% fails R1 with probability 0.93. A setting with a true enrichment-side rate of 6% fails R2 with probability 0.80, and a general excess to 3.5% fails it with probability 0.995.
- **Reported, not gated:** power at linked 0.03 and 0.1; no strata; no clumps; O/E (N2-ev-u, N0-ev-u); the P1 nulls; d = 500 and 2,000 bp.

### G3b: false discovery rate

- **Simulation.** It uses the G3a datasets, not new ones. Per setting, 1,000 random mixes of 15 null and 5 linked (0.03) datasets, with BH at q ≤ 0.05 (mix seed 90,500,000).
  - Pass if the mean share of null pairs among those called is ≤ 0.05 plus Monte Carlo error (one-sided 95% bound).
- **Mock communities (V4, plan step 2.5).** False cross-genome links before and after the chimera module are reported. They are gated only if the mocks carry ≥ 10 true A–B pairs.

## Part 4. Consequences

- **G3a passes:** the net linked test is the read-level test of the paper and of `linkstat test` v0.2.
- **G3a fails in a setting:** the rates and the cause are reported. Read-level p values are not used for claims in data like that setting; the effect sizes with their intervals carry the claims there.
- **G3b fails:** q values are reported as descriptive, with the measured false share.
- G3 does not stop the project. Whether to continue to P3 is the user's decision.

## Feasibility (development seeds; run 7 = the frozen implementation)

Runs 6 and 7 used the G3a design with fewer datasets: 18 settings × {100 under no linkage, 40 at linked share 0.03, 40 at 0.1}, 3,240 datasets each.
- **Run 6** (seeds 84,000,000+, package `2c52ff6`): normal reference; net_w without the margin rule.
- **Run 7** (seeds 86,000,000+, package `100c199`): the text above.

| G3a rule | Run 7, net_w test | Run 7, same-molecule test | Run 6 (net_w; same molecule) |
|---|---|---|---|
| R1 two-sided size, pooled (max setting) | 3.3% (8%) | 4.3% (10%) | 3.7% (6%); 4.3% (8%) |
| R2 enrichment side at 2.5%, pooled (max setting) | 1.6% (4%) | 2.3% (5%) | 2.8% (5%); 2.6% (6%) |
| R3 largest \|mean\| under no linkage | 0.0006 | 0.0020 | 0.0007; 0.0028 |
| R4 lowest coverage of the truth (net_w interval) | 0.90 | (no truth) | 0.60 |
| Power at linked share 0.03 | 0.53–1.00 | 0.33–1.00 | 0.57–1.00; 0.35–1.00 |
| Power at linked share 0.1 | 0.97–1.00 | 1.00 | 1.00; 0.97–1.00 |

G3b in simulation (mixes of 15 null and 5 linked datasets, BH at 0.05): mean false share 0.017 (net_w) and 0.033 (same molecule), so both pass. Run 6 gave 0.035 and 0.048.

What the runs show:
- **Run 6 found the net_w bias** with A in islands under linkage (coverage 0.60–0.80 at 0.1×). The margin rule fixes it: coverage is ≥ 0.90 in every cell in run 7, and 0.935 pooled at 0.1× with islands.
- **The same-molecule maximum of 10%** (10 of 100; soil, 0.1×, A independent) is chance at this size. With runs 6 and 7 pooled, that setting has 8% (16 of 200). Among 36 setting × test cells of 200 datasets, an exactly calibrated test shows a maximum ≥ 8% in 80% of runs. G3a has 500 datasets per setting.
- **The net_w test is conservative.** Rejections on the depletion side are rare when events are few (O is discrete), and the SD of z is 0.96. The simulated genomes appear slightly under-dispersed against the clump model, probably because A and B genes compete for the same space in crowded genomes.
- **R4 is the rule most at risk.** Coverage at 0.1× with islands is 0.935 pooled, with a mean standardised error of −0.28. If that is the true coverage, G3a fails R4 by chance in about 12% of runs.
- **Run time** is 17 s per dataset on one core (0.1×: 25 s; 5×: 5 s; 20×: 22 s), so G3a takes about 5 h on 16 processes.
