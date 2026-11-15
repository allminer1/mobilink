# linkstat (working name)

Tests whether two gene classes sit on the same long read or contig more often than chance, and measures how strongly
they are linked in a way that does not depend on read length.

**Status: pre-alpha (v0.0.1.dev0).** The interface will change, and the name will change before the first release
(`linkstat` is taken on PyPI by an unrelated tool).

## What it computes

- **O/E against a chance model.** Observed co-occurrences against an exposure null: a read of length L can hold a gene
  of length a at about L − a positions. This is a test. It equals 1 under chance at any read length, but under real
  linkage it falls as reads get longer, so it is not an effect size.
- **Net linked fraction at a window w.** For each gene of class A, the distance to the nearest gene of class B on the
  same read is treated as a censored time and estimated by Kaplan-Meier. It is compared with the same curve from
  random anchor positions. The result is 0 under chance and does not depend on read length.

## Inputs

Two tab-separated tables:

- lengths: `molecule_id`, `length`;
- features: `molecule_id`, `start`, `end`, `class`, plus optional annotation columns (0-based, end exclusive).

## Development

```bash
pip install -e ".[dev]"
pytest
ruff check .
```

## Licence

MIT; see `LICENSE`.
