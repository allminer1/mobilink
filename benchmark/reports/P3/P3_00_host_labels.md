# P3 step 3.0: host labels for our soil reads

2026-11-16. Label coverage only; no linkage statistic was computed. Tables: `tables/P3/P3_00_read2bin.tsv`,
`tables/P3/P3_00_soil_kraken.tsv`. This is for the user's decision on the soil host rule, which amendment 04 left
open ("decided before P3").

## Why it matters

The net linked test places its random anchors within host strata. Reads without a host label fall into one
"unresolved" stratum, where genomes rich in resistance genes and rich in mobile elements are mixed, so association
without linkage is not removed there. In simulation, no strata at all gave a net linked fraction of +0.005 to +0.009
under no linkage.

## What each rule labels

"Carriers" are the reads with an ARG, MRG or mobility-core MGE call (the reads the test uses).

| Rule | Carriers labelled, Ridder (4 soils) | Carriers labelled, Pervomaisk (4 soils) | Quality of the labels |
|---|---|---|---|
| A. Amendment default: bin of the contig the read maps to (MAPQ ≥ 20), one assembly | 9–33% | 6–16% | a bin is a genome (CheckM2-filtered) |
| B. As A, over the soil's three assemblies (site co-assembly, then its own two) | 27–53% | 18–29% | as A; bins of different assemblies must be merged by species cluster |
| C. Kraken2 Standard-16, the public-run rule (species, ≥ 10 hit k-mers) | 66% (soil 279) | 63% (soil 304) | poor: a Kraken "species" group is 41% and 25% one MAG; one MAG's reads spread over many groups |
| D. Kraken2 at genus rank | 66% | 63% | poorer still: 31% and 20% one MAG |

- Kraken2 gives more labels, but on soil they are not host labels. The database holds few soil genomes, so reads are
  assigned to whichever relative is present. In the mock communities the same rule was 99.6–99.9% pure, because those
  species are in the database.
- With rule A, 67–94% of a soil's carriers are unresolved; with rule B, 47–82%.

## Options

1. **Rule B as primary (recommended).** It keeps the meaning of the amendment's default ("the bin of the contig the read
   maps to") and labels up to twice as many reads. Unresolved reads stay one flagged stratum. Sensitivity analyses:
   rule A (site co-assembly only), and rule B with Kraken2 genus labels for the reads B leaves unresolved.
2. **Rule A as written.** Simplest, weakest.
3. **Rule B, then Kraken2 for the rest, as primary.** Most reads labelled, but about half of the labels would be the
   noisy kind.

Whatever the rule, the soil read-level results will carry the flag "share of carriers unresolved", and the contig-level
analysis (step 3.4a), where every binned contig has its host, is the stronger evidence for our soils.
