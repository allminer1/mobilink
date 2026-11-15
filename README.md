# linkstat (working name)

Tests whether two gene classes sit on the same long read or contig more often than chance, and measures how strongly
they are linked in a way that does not depend on read length.

**Status: pre-alpha (v0.1.0.dev0).** The interface will change, and the name will change before the first release
(`linkstat` is taken on PyPI by an unrelated tool).

## What it computes

- **Net linked fraction at a window w (the effect size).** For each gene of class A, the distance to the far end of
  the nearest gene of class B on the same molecule is treated as a censored time and estimated by weighted
  Kaplan–Meier. It is compared with the same curve for random anchor positions.
  - It is 0 under chance.
  - In tests on public ONT runs it changed by at most 0.008 when reads were cut from N50 9–14 kb to 3 kb, while the
    naive share fell 2–5×.
- **Molecule-level co-occurrence against a chance model (the test).** Observed molecules with A and B against a
  chance model.
  - Genes cannot overlap, so a molecule carrying B has less room for A. The default test (N2-ev) counts, on each
    molecule, the room an A gene has in the gaps the B genes leave. Within 1 kb bins of that room (and per stratum)
    it is an exact permutation test, conditional on how many molecules carry each class. Length-only nulls are
    biased toward depletion.
  - Reported alongside: N0-ev, a smooth (parametric) approximation of N2-ev, and the registered length-only null.
  - O/E is a test, not a length-invariant effect size.
- **Power and design** (`linkstat power`), and a **simulator with known truth** (`linkstat simulate`).

**Known limits:**

- p values treat molecules as independent. Reads that cover the same locus are not independent, so p values are
  too small at high coverage. Locus clustering is planned; until then, interpret effect sizes and intervals.
- The chance model treats the genes of one class as placed independently. In real metagenomes they come in clusters
  (operons, islands, cassettes) and differ between hosts. When the class the model places is clustered, O/E can lean
  toward depletion. The results report a clustering index for that class and flag values above 2.

## Inputs

Two tab-separated tables:

- lengths: `molecule_id` (or `read_id`), `length`, optionally a stratum column (sample, taxon bin, ...);
- features: `molecule_id`, `start`, `end`, `class`, optionally `scov` (subject cover, %), 0-based, end exclusive.

`linkstat prep` builds the feature table from DIAMOND blastx hits against CARD, BacMet2 and mobileOG-db.

## Quick start (toy data, under a minute)

```bash
# a read-length distribution to sample from (any lengths table works)
linkstat simulate --read-lengths my_reads.lengths.tsv --n-genomes 40 --total-bases 3e8 --linked 0.3 --out-prefix toy
linkstat test --lengths toy.lengths.tsv.gz --features toy.features.tsv.gz --pairs ARG:MGE --out toy.results.tsv
linkstat power --lengths toy.lengths.tsv.gz --features toy.features.tsv.gz --a ARG --b MGE
```

`toy.truth.tsv` holds the true net linked fraction of the simulated community, to compare with `intact_net` in the
results. Settings and input checksums go to `toy.results.tsv.meta.json`.

## Development

```bash
pip install -e ".[dev]"
pytest
ruff check .
```

## Licence

MIT; see `LICENSE`.
