# P2 step 2.5: chimeras and mock communities (V4)

2026-11-16. Design fixed before any result in D48 (`docs/analysis_plan.md`). Development work under the approved P2
plan; the only registered element is the mock part of G3b, read as written in D48 (that reading is for the user to
confirm). Tables: `tables/P2/P2_05_mocks.tsv` (one row per run), `tables/P2/mocks/eval/` (per run: read classes,
links, reference pairs, breakpoints), `tables/P2/soil_adapters/` (our soils).

## Summary

- **Chimeras are rare but vary more than 100-fold between runs.** Reads that join two species: 0.01–0.06% in the two
  rapid-kit runs, 0.03% and 0.42% with native barcoding (our kit), 1.7% and 5.5% with plain ligation, 1.5% in the
  R9.4.1 gut run.
- **False cross-genome links are rare even without any filter.** Pooled over the 7 runs, 92 of 13,686 same-molecule
  A–B links (0.67%) were false before filtering; within 2 kb, 15 of 10,398 (0.14%). In the run with our kit and the gut
  standard (most resistance genes), none were false.
- **G3b mock part (my reading in D48): passes in all 7 runs.** After the chimera module, the upper 95% bound of the
  false share is ≤ 2.4% everywhere (rule: ≤ 5%), for both windows.
- **The overlap rule fixed in D48 is too aggressive.** It flags 12–38% of the clean reads it can judge and would
  discard 13–43% of the true links. A stricter variant (breakpoint ≥ 500 bp from the read ends, ≥ 2 supporting reads)
  flags 0.6–2.2% of those clean reads, catches 73–100% of the two-species chimeras, removes 75–100% of the false links
  in the runs that had them (33% in one run), and loses 0.7–7.5% of true links.
- **Mid-read adapters are a weak signal in current data**: basecallers trim and split at adapters, so only 1–21% of
  two-species chimeras still carry one (73% in the one untrimmed rapid run). In our 8 soils, 0.06% of reads (0.05–0.12%
  per soil) have an adapter or barcode inside the read, mostly the soil's own barcode.
- **Kraken2 host strata are pure**: 99.6–99.9% of reads in a species stratum come from that species; 0.1–1.1% of
  carrier reads are unresolved.
- **Fold-backs** (both strands of one molecule in one read) are 0.07–1.7% of reads. They create no false links but count
  a locus twice on one read.

## Runs

Each run: about 2 Gbp as random blocks of the archived reads (48 blocks via SRA; 8 larger blocks for the two runs held
only at ENA, where fastq-dump needs about 20 min per block; the small rapid run whole). Annotation as for the public
runs (DIAMOND, P0_02c rules), registered classes (ARG and MRG primary rule; MGE = mobility core, D17). References:
Zymo refseq, annotated with the same pipeline. Reads deleted after processing.

| Run | Standard | Kit (from the read ends) | Lab, device | Identity | Chemistry |
|---|---|---|---|---|---|
| ERR14126482 | D6331 gut | native barcoding (NB13), untrimmed | Nanomotif study, PromethION | 0.992 | R10.4.1 (Dorado sup v5.0.0) |
| SRR33364450 | D6322 HMW | native barcoding (NB24), partly trimmed | A*STAR, GridION | 0.982 | R10.4.1 |
| ERR12983023 | D6322 HMW | ligation, untrimmed | Nanomotif study, PromethION | 0.985 | R10.4.1 |
| SRR34704835 | D6322 HMW | ligation, untrimmed | AIT, MinION | 0.935 | R9-like identity |
| SRR28476579 | DNA standard | rapid, untrimmed | compression study, MinION | 0.936 | R9-like identity |
| SRR33788114 | D6300 cells | rapid, trimmed | extraction benchmark, PromethION | 0.954 | unclear |
| SRR23926926 | D6331 gut | barcoded, trimmed | HKU, GridION | 0.958 | R9.4.1 (study title) |

## Read classes against the reference

Alignment segments with MAPQ ≥ 20 and ≥ 300 bp. "Two species" is certain; the total is that rate divided by
(1 − Σ q²), the chance that two random molecules come from different species.

| Run | Two species | Estimated total | Fold-back | Other non-collinear |
|---|---|---|---|---|
| ERR14126482 (NBD) | 0.42% | 0.49% | 0.48% | 0.88% |
| SRR33364450 (NBD) | 0.03% | 0.04% | 1.12% | 0.21% |
| ERR12983023 (ligation) | 1.70% | 2.42% | 1.65% | 0.77% |
| SRR34704835 (ligation) | 5.55% | 8.20% | 0.12% | 1.17% |
| SRR28476579 (rapid) | 0.06% | 0.06% | 0.10% | 0.43% |
| SRR33788114 (rapid) | 0.01% | 0.02% | 0.07% | 0.19% |
| SRR23926926 (R9.4.1) | 1.50% | 1.72% | 1.56% | 0.83% |

"Other non-collinear" (one species, out of order) is mostly larger than the random-joining rate predicts, so it is
largely repeats and split alignments, not chimeras.

## Links: precision before and after the chimera module

Read-level A–B links for ARG:MGE, ARG:MRG and MRG:MGE. True: both calls in one reference segment. False: calls in
segments that are not collinear. Unknown links (a call outside every segment) are 0–1% and left out.

| Run | 2 kb: false / known | Same molecule: false / known | Upper 95% after the D48 module (2 kb; molecule) | G3b-V4 (D48 reading) |
|---|---|---|---|---|
| ERR14126482 | 0 / 566 | 0 / 913 | 0.006; 0.004 | pass |
| SRR33364450 | 0 / 1,771 | 2 / 2,041 | 0.002; 0.002 | pass |
| ERR12983023 | 8 / 2,759 | 52 / 3,517 | 0.004; 0.003 | pass |
| SRR34704835 | 6 / 1,519 | 23 / 2,209 | 0.004; 0.003 | pass |
| SRR28476579 | 0 / 144 | 1 / 159 | 0.024; 0.021 | pass |
| SRR33788114 | 0 / 3,233 | 9 / 4,366 | 0.001; 0.003 | pass |
| SRR23926926 | 1 / 406 | 5 / 481 | 0.013; 0.011 | pass |

Every run spans 68–76 reference A–B pairs within 2 kb (≥ 10, so all are gated in the D48 reading).

## The chimera module

Reference-free flags, scored against the read classes:
- **Adapter inside the read** (P2_05e): false-flag rate 0.07–0.20%, specific. Sensitivity to two-species chimeras
  1–21%, except 73% in the untrimmed rapid run.
- **Overlap breakpoint** (carrier ava overlaps): a cluster of overlap ends inside both reads that no overlap spans.

Rates are among reads that have overlaps (77–97% of carriers).

| Overlap rule | Clean reads flagged | Two-species chimeras flagged | True links lost (molecule; module incl. adapter) |
|---|---|---|---|
| D48 (≥ 200 bp from ends, 1 supporting read) | 12–38% | 82–100% | 13–43% |
| ≥ 500 bp from ends, ≥ 2 supporting reads | 0.6–2.2% | 73–100% | 1.2–7.5% |
| ≥ 1,000 bp from ends, ≥ 1 supporting read | 0.9% (rapid run) | 64% (rapid run) | (computed for one run only) |

Most false breakpoints sit 200–400 bp from a read end, where ONT read ends are frayed. Reads with no overlap cannot
be judged (3–23% of carriers).

## Our soils (mid-read barcodes, D38)

All 8 soils were sequenced on two runs, with barcode NB01–NB08 on run 2780af39 and NB09–NB16 on run cee025ed; both
barcodes of a soil count as its own (`tables/P2/soil_adapters/soil_barcodes_all.tsv`). Full scan of every read for
NB01–NB24 with flanks, the ligation adapter core and the rapid primer (`tables/P2/P2_05_soil_adapters.tsv`).

| Soil | Reads | Any hit inside the read | Own barcode | Another barcode | Ligation adapter | Carriers with a hit inside |
|---|---|---|---|---|---|---|
| 278 | 2,436,725 | 0.056% | 0.045% | 0.004% | 0.010% | 0.048% |
| 279 | 1,022,658 | 0.045% | 0.030% | 0.006% | 0.012% | 0.056% |
| 280 | 1,470,590 | 0.048% | 0.036% | 0.003% | 0.012% | 0.078% |
| 281 | 3,334,828 | 0.063% | 0.048% | 0.005% | 0.013% | 0.066% |
| 304 | 3,718,501 | 0.046% | 0.037% | 0.003% | 0.008% | 0.059% |
| 307 | 5,329,993 | 0.064% | 0.056% | 0.003% | 0.009% | 0.059% |
| 308 | 2,624,845 | 0.069% | 0.061% | 0.004% | 0.007% | 0.115% |
| 310 | 1,675,677 | 0.118% | 0.102% | 0.005% | 0.015% | 0.170% |

- Pooled: 0.063% of reads and 0.078% of carrier reads carry an adapter or barcode inside the read. Dorado had already
  split 2.8% of reads at internal adapters (D38), so these are what read splitting missed.
- Most hits are the soil's own barcode: two molecules of the same sample read as one. Another sample's barcode (a
  cross-sample chimera) is 0.003–0.006% of reads.
- So mid-read adapters can flag only a tiny share of our reads. The mocks put the real chimera rate of native-barcoding
  runs at 0.03–0.42%, so most chimeras in our data would carry no adapter; the overlap flag is the one that matters.

## Reading

- For ONT data from current pipelines, chimeras do not threaten the link results: the run with the most chimeras
  (5.5% of reads) had 1.0% false same-molecule links, and the highest false share was 1.5% (52 of 3,517). False links need a chimera that joins an A fragment to a B fragment, which
  is rare next to the true links in the same data.
- Kit is not the whole story: the two ligation runs differ 3-fold, and the two native-barcoding runs 14-fold. The
  rate is a property of the run (library and run conditions), so a tool should measure it per dataset rather than
  assume it per kit.
- The module as fixed in D48 costs far more true links than it removes false ones. The stricter variant is the sensible
  default for v0.2 (development choice, below).
- Read accuracy drives recall: both genes of a spanned reference pair were called on 92% of reads at 99.2% identity
  and 56% at 93.6% identity.

## Decisions and next steps

- **For the user:** confirm (or change) the D48 reading of the frozen mock rule. All 7 runs pass it either way: the
  false share stays below 2.5% at its upper bound even before any filter.
- **For the user, before P3:** the registered test has no chimera filter. I propose the primary P3 analysis without a
  filter (as registered) and the filtered analysis as a sensitivity analysis. Applying the filter in the primary
  analysis would change a registered definition.
- **Delegated, for v0.2:** a `linkstat chimeras` command that writes the flag table (adapter inside the read; overlap
  breakpoint ≥ 500 bp from the ends with ≥ 2 supporting reads; fold-back by self-alignment), for use with
  `linkstat test --exclude`; the per-run chimera audit in the output.
