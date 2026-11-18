# Amendment 03: excluded-volume null and gate G2

Status: **approved by the user on 2026-11-15 and frozen the same day, before any G2 statistic.**
- Part 1 (D30): "Adopt corrected null".
- Part 2 (D31): "Approve G2 v2". This version replaces draft v1, which the user approved (D29) but which could not be run as written (section 0).
- **Frozen code:**
  - the package at commit `7b665ca7bc62e8553b927556367f6fdbffa0128d` (`frozen/amendment03/PACKAGE_COMMIT`);
  - the G2 scripts in `frozen/amendment03/` (sha256 in `frozen/amendment03/SHA256SUMS`), run by `scripts/P1_15_g2_driver.sh`, which checks both.
- The sha256 of this file is in `docs/preregistration_amendment_03.sha256`.
- The registered pre-registration (`docs/preregistration_H1-H4.md`) and amendments 01 and 02 are unchanged; this amendment changes the registered read-level null (Part 1).

## 0. Why v1 changed before any G2 run

Two things surfaced while setting up G2b. No G2 statistic was computed and no G1b run was analysed with these tests.
Evidence: `reports/P1/step1.15_feasibility.md`, `tables/P1_15_*`.

1. **The density setting "rare = ARG vs MRG on our soils" is not feasible.** It gives about 0.001 expected
   co-occurrences per dataset; a test there never rejects, so there is nothing to calibrate.
2. **Excluded volume.** Genes cannot overlap, so a read carrying B genes has less room for A. Every length-based null
   ignores this: the registered analytic null N0 and the length-binned permutation null N2. In the genome simulator
   with no linkage they give O/E 0.72–0.85 and reject 12–91% at a nominal 5%.
   - On real runs (P0 public runs, development data), the correction lowers the ARG vs mobility-core expectation by 13–28%.
   - Sewage run SRR32022307 flips from "depleted" (O/E 0.95) to "enriched" (1.13).
   - The net linked fraction is not affected: its random anchors already avoid B.

## Part 1. Read-level null with the excluded-volume correction (changes a registered definition)

- **Registered:** "Primary test: O/E under the analytic exposure null (length only)."
- **New primary read-level null, N0-ev.**
  - A's exposure on each molecule counts only the start positions where an A call fits in the gaps the B calls leave:
    - (g − a)+ in a gap between two B calls;
    - (g − f·a)+ in a gap next to a molecule end;
    - (L − (2f − 1)·a)+ on a molecule without B;
    - averaged over the A call-length distribution, with f = 0.8.
  - The null conditions on where B lies and on the number of molecules with A, per stratum. The test is the conditional double saddlepoint (D28).
  - Code: `linkstat.nulls.analytic.excluded_volume_exposure`, `cooccurrence_test(exposure_a=...)`.
- **The same correction applies to N1** (per stratum) **and to N2:** permutation within 1 kb bins of A's exposure instead of molecule length.
- **The registered length-only null is still computed and reported** next to it, so the change is visible.
- The correction is not symmetric in A and B. The class whose positions are conditioned on is the one carried by more molecules in each dataset; for example, MGE for ARG:MGE, and ARG for ARG:MRG in our soils. The co-occurrence count O is the same either way.

## Part 2. Gate G2

### G2a: analytic expectation against the permutation expectation (real data)

- **Data:** the 8 G1b runs (amendment 02 Part C), native reads, primary rule; pairs ARG:MGE (mobility core) and ARG:MRG.
- **Statistic:** ratio = E(N0-ev) / E(N2-ev), per run and pair. Cases with E(N2-ev) ≥ 10 count.
- **Pass:** at least 75% of the counting cases have |ratio − 1| ≤ 0.05, **and** the median |ratio − 1| ≤ 0.05.
- **If it fails:** N2-ev, which is exact under its own permutation null, becomes the primary read-level test; N0-ev is reported as its smooth approximation.

### G2b: type I error under no linkage (simulation)

- **Generator:** `linkstat.simulate`; no linkage; B in islands at the median mobility-core density of the G1b runs (5 islands + 12 single genes per genome, about 32 per 4 Mb); whole and 80% partial calls.
  - Depth is 0.1 × the community size, which puts about 4–10% of covered loci on two or more reads.
  - Strata: the source genome of each read, standing in for the host strata N1 will use (taxon bins).
- **Settings, 3 × 2 = 6:**
  - Read lengths from soil 279 (N50 about 5 kb), wastewater SRR32022521 (10.1 kb) and gut SRR23278673 (14.0 kb).
  - A density 24 per 4 Mb (wastewater SRR32022521), 200 genomes; or 6 per 4 Mb (G1b range), 800 genomes.
  - Expected co-occurrences, measured on 2–3 datasets per setting before freezing: 18–28 with the wastewater and gut read lengths, and 6–8 with the soil read lengths, whose short reads carry fewer genes. Draft v2 said "about 15–20" for all settings; the rules above are unchanged. At E ≈ 6–8 the 3% lower limit may bind, because exact-type tests of small counts are conservative.
- **Datasets:** 4,000 per setting, seeds 20261116 + setting × 100,000 + dataset; only datasets with E ≥ 5 count. With 4,000 datasets a correct test leaves 4–6% by chance in about 0.4% of comparisons.
- **Pass:** in every setting, the two-sided rejection rate of the conditional p value at α = 0.05 is within 3–6%, and the one-sided (enrichment) rate is ≤ 6%.
  - The lower limit of 3% allows for the conservatism of exact-type tests of small counts.
  - The mid-p rates are reported.
- **Reported, not gated:**
  - depth 0.5× and 2× (pseudo-replication; P2's locus clustering);
  - without genome strata (host heterogeneity);
  - the registered length-only nulls N0 and N2;
  - the exposure-model generator of step 1.11.

## Part 3. Consequences

- **G2 does not stop the project** (project plan section 8).
- **G2 passes:** N0-ev (conditional test; mid-p reported for small counts) is the read-level test of the paper and of `linkstat test` v0.1, and the registered null is reported alongside.
- **G2b fails in a setting:** the rates and the cause are reported, and the paper restricts the analytic test to the settings that passed. N2-ev is used where it is exact.
- **Read-level p values stay invalid at high coverage** until P2 adds locus clustering (unchanged from the handoff).

## Part 4. Timeline

G2 runs right after the freeze and takes about 1–2 hours at 16 threads; the G2 memo follows.
