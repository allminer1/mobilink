# Method

This note describes what `mobilink test` computes. The full derivation, the simulations and the pre-registered
validation are in the paper and in `benchmark/`.

## The problem

Take class A (for example ARG) and class B (for example MGE). The usual statistic is the naive share: the share of
A-carrying reads that also carry a B. A longer read has more room for a B, so the naive share grows with read length
even when A and B are placed independently. In simulations with no linkage it rose from 0.07 to 0.38 as the read N50
went from 3 to 24 kb.

## The measure: net linked fraction

1. **Units.** A calls within 1 kb of each other with no B between them form one A unit, so an operon counts once.
2. **Distance with censoring.** For each side of each unit, take the distance to the far end of the nearest B on
   that side. If the read ends first, the observation is censored at the point where a whole B could no longer fit.
3. **Kaplan–Meier.** F_A(w) is the estimated share of sides with a whole B within *w*. Censoring makes it independent
   of read length: a short read contributes only up to where it ends.
4. **Chance.** The same estimate for random anchors gives F_R(w). An anchor has the length of an A call and a random
   position on a read of the same host, and it may not overlap a B call.
5. **Net linked fraction.**

   net_w = 1 − (1 − F_A(w)) / (1 − F_R(w))

   It is 0 under chance, and it equals the share of A units that have a B within *w* because of linkage and not by
   chance. The default *w* is 2 kb.

`net_mol` is built the same way for "a B anywhere on the same molecule". It depends on molecule length and is meant
for contigs and for tests, not for comparison between datasets.

Units that end within 1 kb of a molecule end are left out (and anchors likewise), because a unit cut by the read end
looks complete but is not.

## Hosts

Genomes rich in A are often rich in B. This makes A and B co-occur on reads without being linked on the genome.
mobilink therefore places anchors within the host of each read: Kraken2 species for reads, the bin for contigs, or
any label you give. Each host gets anchors in proportion to its A units.

## The test

The null hypothesis is "no linkage": given where the A units are, B is as near as it is to random positions of the
same host.

- **Statistic.** The number of sides with a B within *w* (or of units with a B on the molecule), minus its
  expectation from the anchors of the same host.
- **Clumps.** Reads that cover the same genomic locus see the same genes, so they are one observation. mobilink
  joins A units into loci from the all-against-all overlaps of the carrier reads, and computes the variance of the
  statistic under the null with whole clumps as units. For contigs, all A units of a contig are one clump.
- **p value.** From the null distribution of the clump model by Monte Carlo (9,999 draws, mid-p), two-sided.
- **Several pairs.** Benjamini–Hochberg q values over the pairs of a run.

Intervals of the effect sizes come from a delete-one-host jackknife.

## What was validated

| Check | Result |
|---|---|
| False positives under no linkage (9,000 simulated datasets, nominal 5%) | 5.0% (`net_w` test) and 4.4% (same-molecule test) |
| Coverage of the true value by the 95% interval (7,200 simulated datasets with linkage, plus the null ones) | at least 0.914 in every setting |
| Simulated truth 0.125, read N50 18.5, 9 and 4 kb | 0.131, 0.130, 0.135 |
| 14 public runs, reads cut to 3 kb | naive share 1.8- to 7.7-fold lower; net linked fraction changed by ≤ 0.006 |
| Contigs, no linkage, genes in clusters or uneven density | label permutation rejects in 67–100% of datasets; mobilink with contig clumps in 0–3% |

## Known behaviour

- After cutting reads with `mobilink oe --fragment`, the estimate is slightly lower (0.5–6% in relative terms, in all
  8 public runs with a net linked fraction of at least 0.06). The cut keeps a gene when 80% of it remains on a piece,
  but the weights assume whole genes, so long genes are over-represented, and long ARG are less often linked to MGE.
  With whole-gene retention the change was centred on zero (−1.3% to +2.7%). `mobilink test` leaves out units near
  molecule ends, so its retention matches its weights.
- With fewer than 30 sides reaching *w*, no estimate is given.
- With no observed event and less than one expected, no estimate is given and the row is flagged "no events".
