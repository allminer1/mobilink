# Changelog

## 0.5.0 (unreleased)

- The package is renamed from the working name `linkstat` to `mobilink` (command, import name, and the version key
  in the `.meta.json` files).
- `mobilink test --molecule-clumps` (API: `loci="molecule"`): all A units of one molecule are one clump. Use it for
  contigs.
- `mobilink test --loci`: only overlaps among molecules with a class-A call are loaded. Results do not change; a
  large overlap file needs much less memory.
- `mobilink power --bases`: the help states the unit (bases).
- Documentation: worked example, real-data guide, output columns, method note; helper tables in `resources/`.
- Benchmark layout (`benchmark/`), CI workflow, Bioconda recipe draft.

## 0.3.0.dev0 (commit 9819182)

The version used for the pre-registered analyses (amendments 05 and 06). Same-molecule ratio with a 90% interval,
pooled runs, filtered overlap reader.

## 0.2.0.dev0 (commit cf0db83)

`test` runs the registered net linked test (amendment 04); chimera flags; reporting fixes.

## 0.1

Molecule-level O/E tables and the first net linked estimator.
