# G1b memo: does the fixed net linked fraction stay stable when reads are shortened?

Date: 2026-11-15. Gate G1b from amendment 02 (Part A approved by the user on 2026-11-14; Part B fixed 2026-11-14;
Part C chosen by the user on 2026-11-15). Every rule below was fixed before any G1b statistic was computed.

## 1. Verdict

**G1b PASSES under the pre-registered rule.**
- All 8 selected runs count and all 8 meet G1b, from 3 environments: pig manure 3, human gut 3, wastewater 2.
- The rule needed at least 6 counting runs from at least 3 environments, with at least two-thirds of them meeting G1b and the meeting runs spanning at least 3 environments.
- Frozen gate output: `tables/P1_08_g1b_gate.tsv`.

Under D16, a pass means P1 core statistics and P2 continue. A4 leaves the next step to the user (section 7).

## 2. What was tested

- **Class pair:** ARG:MGE, with the mobility-core MGE class (D17), all ARGs, primary thresholds, w = 2 kb.
- **Comparison:** within each run, the native reads (subsample N50 ≥ 9 kb) against the same reads cut in silico to N50 3 kb.
  - The cuts use the P0 breakpoint rates, averaged over 20 breakpoint draws.
  - Intervals: read bootstrap with 500 replicates, widened at 3 kb by the between-draw sd.
- **Estimator:** the frozen Part B code (`frozen/partB/`; checksums verified by the driver before it ran).
- **A run counts if:**
  - its subsample N50 is ≥ 9 kb;
  - at least 30 sides reach w at both lengths;
  - the half-width of the 3 kb interval is ≤ 0.10.
- **A run meets G1b if:**
  - the naive share ratio, native / 3 kb, is ≥ 1.5;
  - |native − 3 kb| ≤ max(0.02, 0.25·|native|);
  - the two 95% intervals overlap.

## 3. Runs (amendment 02 A2, A2.1 and Part C)

- **Screen:** 240 eligible runs in 87 independent units. All 240 were probed for read length (3 blocks × 1,000 reads); 104 needed a retry after a time-out (D21).
- **Selection walk (Part C, `public/screen/g1b_selection_walk.tsv`):**
  - Candidates are walked in probe-N50 order, at most 3 per project where an environment has only one project, and one run per patient, site or sample.
  - A run is accepted only if its 3 Gbp subsample has N50 ≥ 9 kb.
- **Replaced runs** (subsample N50 < 9 kb; no statistic computed on them):
  - pig SRR34068984 (5,275), SRR34069004 (7,499), SRR34068991 (7,941);
  - wastewater SRR32022305 (8,828), SRR32022409 (8,862), SRR32022312 (8,739).
- **Dropped environment:** chicken caecum; no run reached 9 kb.
- **Why the walk was needed:** the probe overstated N50 for the runs ranked top. It was off by up to 9 kb in the MinION pig project, where a few very long reads inflate a 3,000-read sample. Wastewater and gut subsamples were within 0.5 kb of their probes.
- **Block checks:** every one of the 48 blocks in each run arrived complete (`public/<run>/block_check.tsv`).

## 4. Results per run (mobility core; `tables/P1_07_g1b_MC.tsv`)

| Run | Environment, unit | N50 | ARG reads | Sides reaching w (native / 3 kb) | Net linked native (95% CI) | 3 kb (95% CI) | 2 kb | Δ 3 kb / tolerance | Rel. change 3 kb | Naive ratio (native → 3 kb share) | Counts / meets |
|---|---|---|---|---|---|---|---|---|---|---|---|
| SRR34068997 | pig manure, sg_17 | 9,149 | 603 | 788 / 315 | 0.060 (0.042–0.075) | 0.059 (0.028–0.089) | 0.058 | 0.0009 / 0.020 | −1.5% | 2.59 (0.325 → 0.125) | yes / yes |
| SRR34068996 | pig manure, sg_18 | 11,087 | 268 | 289 / 90 | 0.099 (0.071–0.132) | 0.107 (0.041–0.174) | 0.098 | 0.0082 / 0.025 | +8.3% | 2.44 (0.321 → 0.131) | yes / yes |
| SRR34068990 | pig manure, sg_23 | 9,130 | 189 | 272 / 81 | 0.058 (0.030–0.082) | 0.058 (−0.005–0.121) | 0.053 | 0.0002 / 0.020 | +0.4% | 2.89 (0.328 → 0.113) | yes / yes |
| SRR23278673 | human gut, patient 23 | 14,041 | 1,014 | 1,477 / 394 | 0.152 (0.137–0.167) | 0.146 (0.115–0.176) | 0.153 | 0.0063 / 0.038 | −4.1% | 2.62 (0.602 → 0.230) | yes / yes |
| SRR23278681 | human gut, patient 21 | 12,322 | 1,011 | 1,897 / 611 | 0.092 (0.079–0.110) | 0.090 (0.065–0.115) | 0.096 | 0.0018 / 0.023 | −1.9% | 2.02 (0.306 → 0.151) | yes / yes |
| SRR23278721 | human gut, patient 28 | 11,128 | 681 | 1,022 / 318 | 0.050 (0.038–0.063) | 0.051 (0.026–0.077) | 0.048 | 0.0012 / 0.020 | +2.3% | 4.98 (0.357 → 0.072) | yes / yes |
| SRR32022521 | wastewater, site KC | 10,124 | 12,826 | 23,816 / 8,154 | 0.013 (0.010–0.016) | 0.012 (0.006–0.017) | 0.013 | 0.0007 / 0.020 | −5.7% | 2.57 (0.153 → 0.060) | yes / yes |
| SRR32022302 | wastewater, site E | 9,616 | 1,428 | 1,854 / 578 | 0.072 (0.060–0.084) | 0.071 (0.047–0.094) | 0.071 | 0.0019 / 0.020 | −2.6% | 2.08 (0.276 → 0.133) | yes / yes |

**In short:**
- **Naive share:** falls 2.0–5.0× when reads are cut from native to N50 3 kb.
- **Net linked fraction:** moves by at most 0.008 in absolute terms (−5.7% to +8.3%; median −1.7%) at 3 kb, and by at most 0.005 at 2 kb (−8.6% to +4.3%; median −1.4%).
- **Tolerance floor:** all 8 runs would also meet a stricter tolerance of 0.25·|native| with no 0.02 floor. The floor did not decide any run.

## 5. Secondary results (not part of the gate)

- **Registered secondary** (median relative change over counting runs with |native| ≥ 0.10): −4.1% at 3 kb and +0.7% at 2 kb. By the A3 wording rule this allows "approximately length-invariant above N50 3 kb".
  - **Caveat:** under the mobility core, only one run (SRR23278673) has net linked ≥ 0.10, so this secondary rests on a single run. On its own it is too thin to support the wording.
  - Two descriptive figures, not registered, back it up: the median over all 8 runs is −1.7% at 3 kb, and the D3 result below.
- **D3 class** (the registered secondary class; frozen gate run on the D3 table, `tables/P1_08_g1b_gate_D3_secondary.tsv`):
  - PASS: 7 runs count, and all 7 meet.
  - SRR34068996 does not count, because its 3 kb half-width is 0.1035 against a limit of 0.10.
  - Registered secondary: median −0.9% at 3 kb over 6 runs with |native| ≥ 0.10; +1.8% at 2 kb.
  - Largest absolute change: 0.011 at 3 kb and 0.013 at 2 kb.

## 6. How far this goes

- **What it shows:** on new runs, with the estimator frozen beforehand, cutting reads from N50 9–14 kb down to 3 kb (or 2 kb) leaves the net linked fraction unchanged within a few percent. Over the same change, the naive share falls 2–5×. This is the main claim of the paper, and it held in every run.
- **New runs, not new projects.** All 3 projects were also P0 sources (PRJNA929328, PRJNA1207805, PRJNA1279599), although every G1b patient, site and sample is new. The one new project, PRJNA1414918 (pig farms), had no run with N50 ≥ 9 kb. Generalising across studies is still open; it is part of P2.
- **The fragmentation is simulated, not real.** The short versions are the same molecules cut in silico. Real short-read runs also differ in DNA extraction, chemistry and basecalling. That comparison needs paired real data (P2 or P3).
- **Small, noisy runs.** The pig runs have 189–603 ARG reads, so their 3 kb intervals are wide (half-width 0.03–0.07). The interval-overlap condition is therefore easy to meet; the stability condition does the real work.
- **Low absolute values.** Under the mobility core, net linked is low (0.013–0.152). The test shows stability, not a large linkage signal.
- **Read length.** Our own soils (N50 about 5 kb) are below the 9 kb used here. G1b supports cutting down to 3 kb, but it says nothing about native runs shorter than 9 kb compared with each other.
- **Process deviations, all before any statistic:**
  - D21: probe retry.
  - D23: estimator held when the A2 selection could not reach a verdict.
  - D24 / Part C: the user's selection rule.

  Each is logged in `docs/analysis_plan.md` section 9, with hashes in `docs/preregistration_amendment_02.part{A,B,C}.sha256`.

## 7. Decision for the user

Under D16 and A3, G1b passing means continuing with **P1 core statistics** (`docs/analysis_plan.md` section 4) and then **P2**:
- analytic and permutation nulls;
- power;
- locus or cluster bootstrap;
- cross-study tests.

The backup project is not triggered.

Points the user may want to settle as P1 starts:
1. **Wording of the length claim.** The registered rule allows "approximately length-invariant above N50 3 kb". Because the mobility-core secondary rests on one run, I suggest also reporting the all-run and D3 figures and keeping "approximately" until P2 adds real cross-run comparisons.
2. **A cross-study check on real short runs.** P2 should include pairs of runs from the same sample sequenced at different read lengths, if such pairs exist in the archive, alongside the Go & Shrestha 2026 datasets (novelty check).

## 8. Files

- Selection:
  - `public/screen/g1b_selected.tsv` (Part C);
  - `public/screen/g1b_selection_walk.tsv`;
  - `public/screen/g1b_selected_A2.tsv` (the A2 selection, not used).
- Per-run results: `tables/P1_07_g1b_MC.tsv` (primary) and `tables/P1_07_g1b_D3.tsv` (secondary).
- Gate:
  - `tables/P1_08_g1b_gate.tsv`;
  - `tables/P1_08_g1b_gate_D3_secondary.tsv`.
- Logs: `logs/P1_07_driver.log`, `logs/P1_07_g1b_MC.log`, `logs/P1_07_g1b_D3.log`, `logs/P1_08_g1b_gate.log`.
- Amendment: `docs/preregistration_amendment_02.md`, Parts A, B and C, with sha256 files.
