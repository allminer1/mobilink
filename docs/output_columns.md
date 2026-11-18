# Output of `mobilink test`

One row per pair A:B. "Unit" means an A unit: A calls within 1 kb of each other with no B call between them
(`--units-d`). "Side" means one side (left or right) of a unit. "Anchor" means a random position of the same length
as an A call, placed on the molecules of the same host; anchors give the chance model.

## Input summary

| Column | Meaning |
|---|---|
| `class_a`, `class_b` | The pair |
| `molecules` | Molecules analysed |
| `excluded_molecules` | Molecules removed by `--exclude` |
| `A_calls`, `B_calls` | Calls of each class, after ORFs called as both classes were removed |
| `n_dual` | ORFs called as both A and B, removed from this pair |
| `units` | A units |
| `unit_loci` | Clumps of units (units that cover the same genomic locus count once) |
| `loci_conflict_share` | Share of units in clumps whose overlaps disagree on orientation |
| `hosts_with_A` | Hosts that carry at least one A unit |
| `hosts_effective` | Effective number of hosts (lower when a few hosts hold most A units) |
| `a_molecules_unresolved` | Share of A-carrying molecules without a host label |
| `anchors` | Random anchors used for the chance model |
| `w`, `a_units_d`, `margin` | Settings: window, unit distance, and the distance a unit must keep from molecule ends |
| `sides_all`, `sides`, `sides_w` | Sides in total, sides used, and sides whose molecule reaches *w* beyond the unit |
| `call_cover_all` | Share of all bases covered by A or B calls |

## Linked within *w* (`net_w`)

| Column | Meaning |
|---|---|
| `net_w` | Net linked fraction within *w*: share of A units with a whole B within *w*, beyond chance. 0 under chance |
| `net_w_se`, `net_w_lo`, `net_w_hi` | Standard error and 95% interval (delete-one-host jackknife; bootstrap over molecules without hosts) |
| `net_w_O`, `net_w_E` | Observed and expected number of sides with a B within *w* |
| `net_w_V` | Variance of the observed number under chance, with clumps |
| `net_w_z` | (O − E) / √V |
| `net_w_p` | Two-sided p value (Monte Carlo mid-p from the null distribution of the clump model) |
| `net_w_p_enrich` | One-sided p value for linkage above chance |
| `net_w_q` | Benjamini–Hochberg q value over the pairs of this run |
| `net_w_eff` | Effective events, E² / V. Below 10 the row is flagged |

## On the same molecule (`net_mol`)

The same columns with the prefix `net_mol_`, for "a B call anywhere on the same molecule". Use these for contigs.

| Column | Meaning |
|---|---|
| `mol_share_A` | Weighted share of A units with a B on their molecule |
| `mol_share_R` | The same share for the random anchors (chance) |
| `mol_ratio` | `mol_share_A` / `mol_share_R`. 1 under chance |
| `mol_ratio_lo90`, `mol_ratio_hi90` | 90% interval, for equivalence tests (for example "within 0.8–1.25") |
| `mol_log_ratio_se` | Standard error of the log ratio |

`net_mol` depends on molecule length, because longer molecules have more room for a B. Compare `net_w` between
datasets, not `net_mol`.

## Descriptive O/E (`oe_*`)

Unit-level co-occurrence on the same molecule against an exposure model. These columns are descriptive; skip them
with `--no-oe`.

| Column | Meaning |
|---|---|
| `oe_fixed_class` | The class whose positions are kept fixed |
| `oe_n2u_observed`, `oe_n2u_expected`, `oe_n2u` | Observed, expected and their ratio |
| `oe_n2u_lo`, `oe_n2u_hi`, `oe_n2u_p` | 95% interval and two-sided p value |
| `oe_n0u`, `oe_n0u_p` | The same ratio under a simpler chance model (no excluded volume) |
| `oe_clustering_index_modelled` | Clustering of the modelled class along molecules (1 = none) |

O/E is 1 under chance at any read length, but under real linkage it falls as reads get longer. Do not use it as an
effect size.

## Flags

| Flag | What to do |
|---|---|
| `net_w: no events (E = …; no estimate)` | Nothing observed and less than one event expected. The data cannot answer; this is not evidence of depletion |
| `net_w: few effective events (E^2/V … < 10)` | The p value rests on very few events; treat it as descriptive |
| `sides < 30` | Fewer than 30 sides reach *w*; the estimate is not reliable |
| `few effective hosts (G_eff …; intervals wide)` | One or two hosts hold most A units |
| `few units (< 30 hosts carry A units)` | The jackknife interval rests on few hosts |
| `…% of A molecules in the unresolved host stratum` | Those molecules share one stratum; host control is weaker for them |
| `no host strata (anchors over all molecules)` | No `--kraken` or `--strata-col`: genomes rich in A and B create association without linkage |
| `no locus clumps (each unit occurrence is its own clump)` | No `--loci`: p values are too small when reads cover a locus more than once |
| `molecule clumps (all A units of a molecule are one clump)` | Information: `--molecule-clumps` was used |
| `locus conflicts (> 10% of units; sides merged)` | Overlaps disagree on orientation for many units; the test is conservative |
| `dense classes (> … of all bases called)` | A and B cover a large share of the bases; crowding can create association |
| `class empty`, `test not possible (…)` | No calls of a class, or too few molecules or calls |

## Other files

- `OUT.meta.json`: version, time, settings, and SHA-256 checksums of the lengths and feature tables.
- `mobilink prep` writes `PREFIX.<rule>.features.tsv.gz` (columns `molecule_id`, `start`, `end`, `strand`, `class`,
  `gene`, `label`, `evidence`, `intrinsic`, `pident`, `scov`, `bitscore`, `qlen`, `db`) and `PREFIX.summary.tsv`.
- `mobilink chimeras` writes one row per examined read: `read_id`, `overlap_assessable`, `n_breakpoints`,
  `flag_overlap`, `n_adapters`, `flag_adapter` and `flag` (either).
