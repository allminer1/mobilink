# Gate G3: result (amendment 04, frozen 2026-11-16)

**G3a passes all four rules for both tests. G3b passes for both tests.** By amendment 04 Part 4, the net linked test is the read-level test of the paper and of `linkstat test` v0.2. Benjamini–Hochberg q values are used as registered.

## Run

- `scripts/G3a_run.sh` ran from 2026-11-16 02:49 to 05:21 UTC on 32 processes.
- **Checks, before any G3 seed was used:**
  - the sha256 of `docs/preregistration_amendment_04.md`;
  - `frozen/amendment04/SHA256SUMS`: the scripts, the read-length profiles and the driver;
  - the package at commit `100c199fab…`, with a clean source.
- The frozen checksums still verify after the run.
- **Datasets:** 16,200, all complete (500 / 200 / 200 per setting and linked share), no failed chunks. Seeds 90,000,000–90,017,899, none repeated. 17.4 s per dataset.
- **Outputs:**
  - per-dataset results: `tmp/G3a/s*.tsv`;
  - per setting: `tables/G3/G3a_settings.tsv` and `tables/G3/G3b_settings.tsv`;
  - verdicts: `tables/G3/G3a_verdict.txt` and `tables/G3/G3b_verdict.txt`.

## G3a: calibration (9,000 datasets under no linkage, 7,200 under linkage)

| Rule | net_w test | Same-molecule test | Rule |
|---|---|---|---|
| R1 two-sided size: pooled (range over settings) | **5.0%** (3.8–6.2%) | **4.4%** (2.8–5.2%) | pooled 2.5–6%, each ≤ 8% |
| R2 enrichment side at 2.5%: pooled (range) | **2.5%** (1.4–4.2%) | **2.4%** (1.4–3.0%) | pooled ≤ 3%, each ≤ 5% |
| R3 largest \|mean\| under no linkage | **0.0005** | **0.0018** | ≤ 0.01 |
| R4 lowest coverage of the truth (net_w interval) | **0.914** | (no truth) | ≥ 0.90 in every cell |
| Power at linked share 0.03 | 0.56–1.00 | 0.35–1.00 | reported |
| Power at linked share 0.1 | 0.99–1.00 | 0.98–1.00 | reported |

- **R4.** Feasibility flagged R4 as the rule most at risk, in the 0.1× islands cells under linkage. Those cells reached 0.92 at linked share 0.03 and 0.97 at 0.1 for soil, and 0.935–0.97 for gut and wastewater. The lowest cell, 0.914, is soil at 0.1× with islands under no linkage.
- **net_w under linkage** matches the truth for independent A (|difference| ≤ 0.0003).
  - With A in islands it is slightly low, mostly in soil:
    - linked 0.03: −0.0004 to −0.0009, for example soil 0.1× 0.0087 against 0.0096 (about 9% of the effect);
    - linked 0.1: up to −0.0012, for example soil 20× 0.0306 against 0.0317 (about 4%).
  - Elsewhere it is within ±0.0003 at 0.03 and up to −0.0008 at 0.1.
  - The intervals absorb it: coverage is ≥ 0.92 in all these cells.
- **Power at linked share 0.03.**
  - The net_w test: 0.98–1.00 at 0.1× with 2,000 hosts; 0.86–0.93 at 5–20× with independent A; 0.56–0.68 at 5–20× with A in islands.
  - The same-molecule test is less powerful: 0.35–0.54 at 5–20× with A in islands.

## G3b: false discovery rate in simulation

- 1,000 mixes per setting of 15 null and 5 linked (0.03) datasets, with BH at q ≤ 0.05 (mix seed 90,500,000).
- Mean share of null pairs among those called:
  - net_w: 0.031 (MC SE 0.0007), mean power 0.69;
  - same molecule: 0.029 (MC SE 0.0007), mean power 0.51.

  Both pass (rule: ≤ 0.05 plus the one-sided 95% Monte Carlo bound).
- The largest setting values are 0.053 and 0.057. They are reported, not gated.
- Mock communities (V4) are still to come in plan step 2.5.

## Development runs and G3a differ in the net_w size

- **The difference.** The amendment's feasibility section called the net_w test conservative: 3.3% pooled in run 7 and 3.7% in run 6. G3a, with five times as many datasets, gives 5.0%, which is nominal (z SD 0.997, 1.000 and 0.992 at 0.1×, 5× and 20×).
  - The gap lies at 5× and 20×, where a dataset has about 7 effective events. There, z > 1.96 occurred in 3.0% of datasets in runs 6–7 together (71 of 2,400) against 4.35% in G3a (261 of 6,000), about 3.3 SE apart.
  - In run 7 at 20× it was 7 rejections of 600 (1.2%), against 144 of 3,000 in G3a.
- **What was checked:**
  1. **Same code.** Run 7 started at 21:20:43, the time of commit `100c199`; the package was last edited at 21:13. Re-running 12 run 7 datasets with the frozen script reproduced every value exactly.
  2. **Seeds.** No seed is repeated within or between runs. The simulator's derived stream (seed + 5,555,555, used for the truth) collides with no other dataset's seed.
  3. **Independence.** There is no lag-1 correlation of z between neighbouring datasets. The spread of rejection rates across settings is binomial or smaller in all runs. Within G3a, datasets 0–99 and 100–499 agree (4.5% and 4.9% at 20×).
  4. **Simulated data.** O, E, V, effective events and truth have the same distributions in run 7 and G3a (Kolmogorov–Smirnov p 0.06–0.95).
- **Reading.** No code, seed or data difference was found. The development runs were small and are not part of the registered result. G3a is the registered calibration and shows the net_w test at its nominal size.
  - The claim "the net_w test is conservative" in the frozen feasibility section is not supported by G3a. It should not be repeated in the paper.
  - This changes no rule and no result: the R1 floor (2.5%) did not bind, and 5.0% would also have passed the original floor of 4%.
- **Design note for v0.2.** Every dataset uses the test's fixed default seed (`nl.SEED`) for anchors and Monte Carlo draws. This is valid for each dataset, but in simulation the anchor and Monte Carlo noise is common to all datasets. G3 cannot show whether that matters. A per-dataset test seed in simulation studies is the cleaner design, to be checked in v0.2 (descriptive, not a gate).

## The earlier draft rules, for transparency

The pass rules were revised before G3, on development data only (D45). They were revised because an exactly calibrated test would often fail the old ones. G3a shows this:
- the draft rule "3–6% in every setting" would have failed net_w (6.2% in wastewater 20× with islands) and the same-molecule test (2.8% in wastewater 20× with independent A);
- the draft rule "enrichment side ≤ 4% in every setting" would have failed net_w (4.2% in wastewater 5× with independent A).

Under the frozen rules both tests pass.

## Next

- **P2, continuing under the approved P2 plan (D36):**
  - step 2.5: chimeras and mock communities (G3b V4);
  - step 2.6: annotation sensitivity;
  - step 2.8: `linkstat test` v0.2 with the net linked test, Kraken2 strata and `--loci`, plus the per-dataset seed check above;
  - then the P2 report.
- **The user decides** whether to continue to P3 after the P2 report.
