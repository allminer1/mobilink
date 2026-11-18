# Amendment 02 to the pre-registration of H1 to H4: mobility-core MGE class and gate G1b

Status: **Part A approved by the user on 2026-11-14** (A2.1 clarifications added the same day, before any G1b probe). Part A is frozen: sha256 of the text above the "## Part B" heading is recorded in `docs/analysis_plan.md` (D18). Part B is fixed at the end of the P1 development
block, before any G1b data are annotated. Nothing in this amendment changes `docs/preregistration_H1-H4.md`
(sha256 `b1212d82…`) or amendment 01. It adds to them.

Background: G1 failed under its fixed rule (D10; `reports/P0/G1_memo.md`). On 2026-11-14 the user chose option B: fix
the estimator in P1, then a second, pre-registered gate G1b on new data (D16). The user also approved a
data-independent "mobility core" MGE class (D17).

## Part A (for approval now)

### A1. Mobility-core MGE class (MC)

A mobileOG-db beatrix-1.6 entry is in MC if all four rules hold. The rules use only mobileOG's own annotation fields
(`scripts/P1_00_mobility_core.py` → `db/mobility_core.tsv.gz`, 69,026 of 775,257 entries,
sha256 `6cbb08eefe9708428eed1e2b1224a81fba36615393e45ee410f404ae37b15b32`):

1. **Major category** is integration/excision or transfer: DNA that moves within or between cells. Phage is reported as
   its own class. Replication/recombination/repair and stability/transfer/defense are out.
2. **Evidence** is Manual or Homology, not Keyword Search.
3. **Minor category** names none of: regulation, partitioning, competence, replication/recombination/repair. These are
   regulators, plasmid maintenance, natural transformation and housekeeping recombination, all with chromosomal
   homologs.
4. **Gene family** is not in a fixed list of chromosomal housekeeping families. The list covers chromosome-dimer
   recombinases (xerC/D), recombination and repair, replication, segregation and division, Clp proteases and
   two-component response regulators; 42 names, listed in the script. After rules 1–3 it removes xerC/xerD (646
   entries), copR (205), parA (8) and ftsK (2).

Thresholds as registered for MGE: identity ≥ 60%, subject cover ≥ 80%.

Honesty note: copR, recA, recX, clpP and nth were first noticed in our soil results (amendment 01). Under rules 1–3, only
copR is still relevant, and it is on the list as a two-component response regulator family. The other 37 names come
from the functional groups above, not from our hits.

**Role.** MC is the **primary MGE class** for G1b and for every confirmatory test from P1 on. The registered D3 class
and phage are reported next to it as secondary classes. All P0 results stay as they were, on D3.

### A2. G1b data (selected before any annotation)

- **Environments (4):** wastewater, human gut, pig farm/manure, chicken caecum.
- **Eligible runs:** ONT shotgun, native DNA, random selection, ≥ 6 Gbp in the archive, not analysed in P0. Excluded:
  the 8 P0 runs and ERR17169941.
- **Projects** (from `public/screen/ena_project_screen.tsv`):
  - wastewater: PRJNA1207805;
  - human gut: PRJNA929328, ERP188941, PRJEB108131;
  - pig farm/manure: PRJNA1414918 (new), PRJNA1279599;
  - chicken caecum: ERP183232.
  - Each environment includes a project not used in P0 where one exists: pig farm (PRJNA1414918) and human gut
    (ERP188941 and PRJEB108131).
- **Choice:** the N50 of every eligible run is probed as in D10 (3 blocks × 1,000 reads). Per environment, the 3 runs
  with the longest probe N50 and N50 ≥ 9 kb are taken, at most 2 from any one project. Ties are broken by a random
  draw seeded with the accession.
  - An environment with fewer than 2 qualifying runs is dropped, and the drop is recorded.
  - R9 runs, such as PRJNA745436, are reported as their own stratum and do not count.
- **Subsample and annotation:** 3 Gbp per run by random blocks, as in P0 (validated in D14), then the P0 DIAMOND
  pipeline.
- **A2.1 Clarifications (D18, 2026-11-14, after approval and before any G1b probe).** The archive labels did not match
  the wording above, so A2 is read as follows. The list is `public/screen/g1b_eligible.tsv` (`scripts/P1_01_g1b_eligible.py`):
  240 runs, 87 independent units.
  - "Shotgun, random selection" is read as in D10:
    - library source METAGENOMIC;
    - strategy WGS, or OTHER where the project describes untargeted shotgun metagenome sequencing (PRJNA1414918);
    - selection RANDOM or unspecified (D10 accepted PRJNA929328 with "unspecified").
  - Mock, reference and synthetic samples are excluded (Zymo standards, "Reference_batch"), and so are isolates.
  - Projects that turned out to have no eligible run: ERP188941 (its ≥ 6 Gbp runs are isolates) and PRJEB108131 (no run
    ≥ 6 Gbp). Human gut therefore rests on PRJNA929328 alone, so at most 2 gut runs can be chosen.
  - **One run per independent unit:** gut patient (sample alias prefix); wastewater site (any season); pig farm; manure
    sample; chicken. Units already used in P0 are excluded: gut patients 15, 16 and 17; wastewater sites CW, I and TP.

### A3. G1b test

The **unit is one run**. P0 showed that pooling runs mixes unlike samples (D15).

**A run counts if all of:**
- subsample N50 ≥ 9 kb;
- at least 30 sides reach w at both native length and random N50 3 kb (as registered);
- the 95% interval of the 3 kb estimate has a half-width ≤ 0.10, so a run cannot pass on noise alone.

**A run meets G1b if all of:** for ARG:MC, all ARGs, primary thresholds, w = 2 kb, with the Part B estimator:
- **naive share:** native / random N50 3 kb ≥ 1.5;
- **net linked stable:** |native − 3 kb| ≤ max(0.02, 0.25 × |native|), the D10 tolerance unchanged;
- **intervals:** the 95% intervals overlap.

The 3 kb estimate is the mean over 20 breakpoint draws (seeds 20261114 + 0…19). Its interval is the Part B interval at
seed 20261114, widened by the between-draw sd.

**G1b passes if:**
- at least 6 runs count, from at least 3 environments; **and**
- at least two-thirds of the counting runs meet G1b; **and**
- runs that meet G1b come from at least 3 environments.

If fewer than 6 runs count, G1b is **inconclusive** and goes back to the user.

**Secondary, reported but not part of the gate:** the median relative change native → 3 kb over counting runs with
|native| ≥ 0.10.
- If it is within ±10%, the paper may say the estimator is "approximately length-invariant above N50 3 kb".
- Otherwise it says "much less length-dependent than the naive share", with the bias quantified.
- The same is reported at 2 kb.

**Consequence (D16):**
- **G1b passes:** continue with P1 core statistics and P2.
- **G1b fails:** switch to the backup project (simulator and benchmark).
- **G1b is inconclusive:** the user decides.

### A4. Timeline

| When | What |
|---|---|
| P1 weeks 1–2 (from 2026-11-16) | Novelty check (PubMed, bioRxiv). Estimator fixes, each validated on simulation with known truth and on the P0 data as development data: empirical chance term, gene-length correction, short-read censoring check, cluster-robust interval. G1b run screen (N50 probes only). |
| End of week 2 | Part B written and fixed (date, sha256); then G1b fetch and annotation start. |
| Week 3 (project week 5) | G1b result and memo. **The user decides on the next step.** |

## Part B: the G1b estimator (fixed 2026-11-14, before any G1b run was annotated)

Status: fixed. The sha256 of this section (from the "## Part B" heading to the erratum) is in `docs/preregistration_amendment_02.partB.sha256`.

Chosen on simulations with known truth and on the P0 data only. Evidence:
- `reports/P1/step1.2-1.6_estimator.md`;
- `tables/P1_02_sim_scores_partial.tsv`, `tables/P1_05_interval_reps.tsv`;
- `tables/P1_03_dev_p0_partial.tsv`, `tables/P1_06_g1b_dryrun_p0_*.tsv`.

**Code and data, frozen.** The G1b run uses only the copies in `frozen/partB/` (read-only); their sha256 sums are in
`frozen/partB.sha256`:

| File | sha256 |
|---|---|
| P1_06_g1b.py (driver; loads only the frozen copies) | `7ea223c539e4a50f5e7da5c89a386f993c599927bfd497ae9f244ba5ee4549df` |
| P1_03_estimators.py | `0129b455e94aa0099f38e55db319781ade52b77e8760d621a3d7c658716d6b7d` |
| P0_04_gonogo.py (same-ORF rule, breakpoints, N50) | `65ec83c632f91ac0c81338313982c8b4c985190fc6c3fcb5518e7237082c1bcd` |
| P0_02c_features.py (hit filter and one call per ORF) | `2dcb916518d3479e29261b5ab7d7126dce902223788139afb26b1d2b60f89c12` |
| mobility_core.tsv.gz (A1, corrected per the erratum) | `b0bcbcaa8ed4733f6e46abb2ded2a7ae497c57c6d35465e235f1b36a636bae92` |

The frozen driver reproduces the development dry run exactly (`DataFrame.equals`, run SRR34068994).

**B1. Classes.**
- **A:** ARG calls under the primary rule (CARD 4.0.1 protein homolog, identity ≥ 80%, subject cover ≥ 80%), all ARGs.
- **B:** mobileOG-db beatrix-1.6 calls at identity ≥ 60% and subject cover ≥ 80%. One call per ORF over all mobileOG hits, as in P0. A call belongs to B if its subject is in the mobility core.
- **Same-ORF rule** between A and B, as registered.

**B2. Distance and censoring.** For each side of each A call:
- event time = distance to the **far end** of the nearest B call on that side;
- censored at the read end;
- a B call **cut by a read end** censors the side instead of ending it. A call is cut if it ends exactly at a read or fragment end, or within 50 bp of the read end with subject cover < 95%.
- Why: P0's rule (distance to the B start, censored at the read end minus the mean B length) loses long B genes near read ends. That loss grows as reads shorten, and it caused the gut short-read bias.

**B3. Chance term:** the same estimate for 20,000 random anchors.
- Anchor lengths are drawn from the A calls; anchors are placed uniformly on the reads where they fit, seed 20261114.
- Anchors overlapping a B call by ≥ 50% are dropped.
- Anchors are fragmented together with A and B.
- net(w) = 1 − (1 − F_A(w)) / (1 − F_anchor(w)).
- Why: the Poisson term exp(−λ_B·w) is biased when MGEs cluster (−0.02 to −0.04 in simulation, −0.04 to −0.08 in dense public samples).

**B4. Gene-length weights.** Each A side, and each anchor side, is weighted by 1 / P(a gene of its observed length lies whole on a read of the current read set), with P(g) ∝ Σ_r max(0, L_r − g). The Kaplan-Meier estimate is weighted.
- Why: fragmentation and short reads lose long genes. Where linkage depends on gene length, this shifts the estimate by up to +0.026 without weights.

**B5. Window and minimum count.** w = 2 kb (registered). NaN if fewer than 30 A sides reach w (unweighted count).

**B6. Fragmentation.**
- P0 breakpoint rates for N50 3 kb and 2 kb, 20 draws each, breakpoint seeds 20261114 + 0…19.
- A call cut by a breakpoint is kept on the piece holding ≥ 80% of it, truncated (the caller's subject-cover rule).
- The reported 3 kb and 2 kb values are the means over the 20 draws.

**B7. Intervals.**
- Percentile 95% bootstrap over reads, 500 replicates. A sides and anchor sides are resampled by read; the bootstrap seeds are fixed in the driver.
- 3 kb and 2 kb: half-width = √(h² + (1.96·sd)²), where h is the bootstrap half-width at breakpoint seed 20261114 and sd is the between-draw sd.
- Why: on 40 re-sequencings of one simulated community, the read bootstrap covered the truth 88–98% of the time at native length and at 3 kb.
  - G1b compares two versions of the same sample, so re-sequencing variance is the right target.
  - A locus bootstrap (98–100%, 2–3× wider) measures generalisation across loci. It is kept for P2's cross-study tests.

**B8. Gate.** Applied per run exactly as in A3, with these values: counts, naive ratio, stability with the D10 tolerance, overlap, meets. Then the G1b rule across runs (A3).
- The secondary columns are the relative changes native → 3 kb and → 2 kb, for runs with |native| ≥ 0.10.
- D3 is reported with the same code (`--b-class D3`) as a secondary class.

**Validation summary (development only).**
- **Simulation, 5 scenarios × 5 seeds, with caller-like partial calls:** mean error against the genome truth within ±0.007 at every N50 from 2 kb to native (±0.0054 at 3 kb).
  - The P0 prototype was off by −0.02 to −0.04.
- **P0 public runs, development data:** every run with N50 ≥ 9 kb meets the A3 columns under the mobility-core and D3 classes. The largest native vs 3 kb difference is 0.006 (mobility core) and 0.007 (D3).

## Erratum to A1 (2026-11-14, D19; before any G1b annotation)

The rule in A1 is unchanged. The list built from it in the first version was wrong.
- **The bug:** `scripts/P1_00_mobility_core.py` took evidence and minor category from mobileOG subject IDs as DIAMOND reports them. Those IDs are cut at the first space, which removes the evidence field in 114,789 entries and cuts some minor categories. Rule 2 therefore dropped every entry whose evidence field was missing, including 81,317 Homology or Manual entries.
- **The fix:** the script now reads the database's own metadata table (`db/mobileog_beatrix1.6/mobileOG-db-beatrix-1.6-All.csv`), which holds the full categories and evidence. It keys entries by the full header cut at its first space, which is how DIAMOND reports them.
- **Corrected list:** `db/mobility_core.tsv.gz` has **109,325** entries (not 69,026), sha256 `b0bcbcaa8ed4733f6e46abb2ded2a7ae497c57c6d35465e235f1b36a636bae92`.
  - Removed by the name list after rules 1–3: copR 438, xerD 343, xerC 334, parA 9, ftsK 4, uvrA 1.
- The Part A text above, and its sha256, are left as approved. This erratum supersedes the count and sha256 given in A1.

## Part C: G1b run selection after the screen (user decision D24, 2026-11-15; before any G1b statistic)

**Why.** The A2 rule picked 6 runs, the A3 minimum. Two facts appeared before any estimate:
- In every environment, all runs qualifying by probe come from a single project, so the cap of 2 per project limited G1b to 6 runs.
- Pig run SRR34068984 has subsample N50 5,275, against a probe N50 of 14,675. The probe of 3,000 reads is inflated by a few very long reads in the MinION runs. A3 counts a run only with subsample N50 ≥ 9 kb.

Under A2 as written, G1b was therefore INCONCLUSIVE before any estimate. The user chose this rule on 2026-11-15, before any ARG, MRG or MGE statistic had been computed on any G1b run. Only read lengths were known.

**C1. Selection (replaces the cap and the order of checks in A2 "Choice"; A2.1 eligibility and units unchanged).**
1. **Candidates:** eligible runs with a complete probe (3 blocks × 1,000 reads) and probe N50 ≥ 9 kb. Within each environment they are ordered by probe N50, longest first, with ties broken by the random draw seeded with the accession as in A2 (`public/screen/g1b_blockprobe.tsv`, 240 of 240 runs probed).
2. **Project cap:** 3 runs per project when all of an environment's candidates belong to one project; otherwise 2.
3. **Walking down the list,** a candidate is accepted if:
   - its unit (gut patient, wastewater site, pig farm, manure sample, chicken) has no accepted run;
   - its project is below the cap;
   - its 3 Gbp subsample has N50 ≥ 9,000. The subsample is taken by random blocks as in A2 (`public/<run>/subsample_stats.tsv`).
4. A candidate whose subsample N50 is < 9,000 is **replaced**: it is not used in G1b, it does not use up its unit, and the walk goes on. Subsamples are produced in list order, only for runs the walk reaches.
5. The walk stops at 3 accepted runs per environment or when the candidates run out. An environment with fewer than 2 accepted runs is dropped, and the drop is recorded.

**C2. Test.** A3 is unchanged and is applied to the accepted runs only. So are the Part B estimator and the frozen gate. The estimator runs once, after the selection is final. Replaced runs are listed in the G1b memo with their subsample N50, and no statistic is computed on them.

**C3. State when this part was fixed.**
- No G1b statistic had been computed.
- Subsample N50 known for 9 processed runs:
  - pig: SRR34068984 5,275; SRR34068997 9,149; SRR34069004 7,499;
  - gut: SRR23278673 14,041; SRR23278681 12,322; SRR23278721 11,128;
  - wastewater: SRR32022521 10,124; SRR32022302 9,616; SRR32022305 8,828.
- Pig runs SRR34068996, SRR34068991 and SRR34068990, the next pig candidates, had been processed, but their N50 had not been looked at.
- Implementation: `scripts/P1_01g_g1b_select_c.py`.
