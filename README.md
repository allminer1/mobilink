# linkstat (working name)

Tests whether two gene classes sit on the same long read or contig more often than chance, and measures how strongly
they are linked in a way that does not depend on read length.

**Status: pre-alpha (v0.2.0.dev0).** The interface will change, and the name will change before the first release
(`linkstat` is taken on PyPI by an unrelated tool).

## What it computes

`linkstat test` runs the net linked test (pre-registered, amendment 04; calibrated in 16,200 simulated datasets at
5.0% and 4.4% size):

- **Effect sizes.** For each A unit (A genes within 1 kb with no B between them), the distance to the far end of the
  nearest B on the same molecule, as a censored time, estimated by weighted Kaplan–Meier and compared with random
  anchor positions on the same molecules:
  - `net_w`: linked within w (default 2 kb); 0 under chance;
  - `net_mol`: on the same molecule; 0 under chance.
  Anchors are placed within host strata (Kraken2 species, or any column you give), because genomes rich in A are often
  rich in B, which is association without linkage. Intervals come from a delete-one-host jackknife.
- **Tests.** Score tests whose variance is computed under the null, with reads that cover the same locus as one clump
  (from minimap2 all-vs-all overlaps among the gene-carrying reads), and p values from the exact null distribution by
  Monte Carlo (mid-p). Benjamini–Hochberg q values across pairs, per test.
- **Flags** for few effective events or hosts, dense classes, missing strata or clumps, and unresolved hosts.
- **Descriptive:** unit-level O/E against an exposure null (columns `oe_*`); the older molecule-level table is
  `linkstat oe`.

Also: `linkstat chimeras` (reference-free flags for chimeric reads: overlap breakpoints, adapters inside reads),
`linkstat power` and a simulator with known truth (`linkstat simulate`).

**Known limits:**

- With one dominant host (for example a sample that is mostly E. coli), there are few independent hosts and loci, and
  the intervals are honestly wide.
- The effect size depends on what counts as A and B (identity thresholds, which mobile-element categories): report the
  class definitions with every estimate.

## Inputs

Two tab-separated tables:

- lengths: `molecule_id` (or `read_id`), `length`, optionally a stratum column (sample, taxon bin, ...);
- features: `molecule_id`, `start`, `end`, `class`, optionally `scov` (subject cover, %), 0-based, end exclusive.

`linkstat prep` builds the feature table from DIAMOND blastx hits against CARD, BacMet2 and mobileOG-db.

## Quick start (toy data, under a minute)

```bash
# a read-length distribution to sample from (any lengths table works)
linkstat simulate --read-lengths my_reads.lengths.tsv --n-genomes 40 --total-bases 3e8 --linked 0.3 --out-prefix toy
linkstat test --lengths toy.lengths.tsv.gz --features toy.features.tsv.gz --pairs ARG:MGE --strata-col genome \
  --out toy.results.tsv
linkstat power --lengths toy.lengths.tsv.gz --features toy.features.tsv.gz --a ARG --b MGE
```

`toy.truth.tsv` holds the true net linked fraction of the simulated community, to compare with `net_w` in the
results. Settings and input checksums go to `toy.results.tsv.meta.json`.

Real data, with Kraken2 labels for all reads (confidence 0) and overlaps among the gene-carrying reads:

```bash
minimap2 -x ava-ont carriers.fa carriers.fa > carriers.ava.paf
linkstat chimeras --paf carriers.ava.paf --out chimeras.tsv          # optional
linkstat test --lengths reads.lengths.tsv --features reads.features.tsv --pairs ARG:MGE --pairs ARG:MRG \
  --kraken reads.kraken.tsv --taxonomy k2db/ktaxonomy.tsv --loci carriers.ava.paf --out results.tsv
```

## Development

```bash
pip install -e ".[dev]"
pytest
ruff check .
```

## Licence

MIT; see `LICENSE`.
