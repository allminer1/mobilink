# linkstat: pre-registered hypotheses and test settings

Status: **fixed on 2026-11-14**, under the user's instruction "i trust you, take decision by you and continue". It was fixed before any linkage result existed: the P0 runs on our own soils had not finished. Its sha256 is recorded in `docs/analysis_plan.md`, section 9. The file is not edited after that; later changes go into a dated amendment section at the end.

Scope: all analyses of public data, and the cross-study claims. The Ridder and Pervomaisk go/no-go runs (P0 steps 0.1 to 0.5) are pipeline checks and are not covered.

## Hypotheses

- **H1.** Across studies, the naive co-localization share (reads with A and B / reads with A) rises with read N50.
- **H2.** Where there is no linkage, O/E stays at 1 at any N50. Where there is linkage, the net linked fraction at a fixed window does not depend on N50.
- **H3.** ARG-MGE linkage is enriched (O/E above 1) in wastewater and clinical samples more than in soil.
- **H4.** ARG-MRG co-carriage is at or below chance in most soils, polluted ones included. This test is two-sided.

## Definitions

| Item | Setting |
|---|---|
| Caller | DIAMOND 2.2.8 blastx, `-F 15 --range-culling --top 10`, default sensitivity, e ≤ 1e-3; uncompressed FASTA input |
| One call per ORF | Within each read and database, hits are taken by bitscore (ties: e-value, then subject ID). A hit overlapping a kept hit by more than 50% of the shorter one is dropped, and the kept hit labels the ORF |
| ARG | CARD 4.0.1 protein homolog model; identity ≥ 80%, subject cover ≥ 80% |
| MRG | BacMet2 experimental entries that name a metal or metalloid (keyed on ID + accession); identity ≥ 80%, subject cover ≥ 80% |
| MGE | mobileOG-db beatrix-1.6 All; categories integration/excision, transfer, replication/recombination/repair, phage; **identity ≥ 60%**, subject cover ≥ 80% |
| Same-ORF rule | A and B calls overlapping by ≥ 50% of the shorter are removed from that pair's test |
| Primary pairs | ARG:MGE, ARG:MRG, MRG:MGE (class level) |
| Primary effect size | Net linked fraction at **w = 2 kb**, reported only where at least 30 sides reach w |
| Primary test | O/E under the analytic exposure null (length only). A cluster-robust interval is added once P2 locus clustering exists; until then, effect sizes are interpreted, not p values |
| Multiplicity | Benjamini-Hochberg across pairs within each level |
| "At chance" | TOST equivalence on log(O/E) with bounds **±log(1.25)** |
| Sensitivity analyses, fixed now | ARG and MRG at 60/80; MGE at 80/80 and 40/80; ARG without the intrinsic templates (vanR_in_vanO_cl, rpoB2, HelR); MGE with mobileOG evidence Manual + Homology only; MGE per category; windows 1 and 5 kb |
| Sample selection | Shotgun native DNA, ≥ 1 Gbp, read N50 ≥ 2 kb, kit and basecaller known, R9 as its own stratum; every exclusion recorded |
| Subsampling | A random subsample of 2 to 5 Gbp per sample (not the first reads of a file), with the seed recorded |

## Reasons for the four choices fixed on 2026-11-14

1. **MGE identity 60%.**
   - MGE is a functional class (transposases, integrases, conjugation and phage proteins) that is very diverse in soil. On the 1 Gbp pilot, 80% identity kept only 1,085 MGE reads per Gbp, against 8,753 at 60%.
   - 40% sits at the floor of DIAMOND's default sensitivity, and it raises the share of the least curated mobileOG entries (Keyword Search) from 23% at 80% to 32% at 60% and 37% at 40%.
   - 60% is the middle ground. The other two thresholds are fixed sensitivity analyses.
2. **Window 2 kb.** Our reads have an N50 of about 5 kb and the in-silico cuts go down to 2 to 3 kb. At 5 kb, the Kaplan-Meier estimate is undefined for most cut levels and short-read datasets, because fewer than 30 sides reach the window. At 1 kb it captures only adjacent genes.
3. **TOST bounds ±log(1.25).** This is the conventional equivalence margin (0.80 to 1.25), and it is narrow enough to call a result "at chance" only when the interval is tight.
4. **Where registered.** In the project repository now, with the date and sha256. Deposit on OSF or Zenodo when the preprint is posted; that is the user's step, because it needs their account.

## Amendments

(none)
