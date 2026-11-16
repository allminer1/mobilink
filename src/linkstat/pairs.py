"""`linkstat test` (v0.2): the registered net linked test for a set of class pairs, one row per pair.

Registered analysis (amendment 04, frozen 2026-11-16; gate G3 passed, D47):
- A is measured as units (calls within units_d = 1 kb with no B call between them); the same-ORF rule removes A and
  B calls that overlap by >= 50% of the shorter (pre-registration H1-H4).
- Effect sizes: net_w (w = 2 kb) and net_mol (same molecule) against random anchors placed within host strata, with
  delete-one-host jackknife intervals (t on G_eff - 1); without strata, a bootstrap over molecules (flagged).
- Tests: null-variance clump score tests with Monte Carlo mid-p (9,999 draws). Clumps are unit loci from minimap2
  ava-ont overlaps among the carrier reads (`--loci`); without them each unit occurrence is its own clump (flagged:
  anti-conservative when reads cover a locus more than once).
- Benjamini-Hochberg q values across the pairs, separately for each test (pre-registration: "across pairs within each
  level").
- Host strata: Kraken2 per-read labels (linkstat.hosts; species, coarser taxa kept, else "unresolved") or any column
  of the lengths table. The share of A-carrying molecules in the unresolved stratum is reported and flagged.
- Chimeras: molecules listed in `exclude` (e.g. the output of a chimera filter) are removed before anything else; the
  count is reported.
Descriptive, not tests: unit-level O/E under N2-ev-u (exact, binned on unit room) and N0-ev-u, with the fixed class of
amendment 03 (the class whose units sit on more molecules).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .features import drop_same_orf
from .hosts import UNRESOLVED
from .nulls.units import cooccurrence_test_units
from .stats import netlinked as nl
from .stats.multiple import bh
from .stats.nettest import loci_from_paf, net_linked_test


def net_test_pairs(
    lengths: pd.DataFrame,
    features: pd.DataFrame,
    pairs,
    strata=None,
    paf: pd.DataFrame | None = None,
    loci=None,
    exclude=None,
    w: int = 2000,
    units_d: int | None = 1000,
    seed: int = nl.SEED,
    n_mc: int = 9_999,
    oe: bool = True,
    same_orf: float | None = 0.5,
) -> pd.DataFrame:
    """lengths: molecule_id, length; features: molecule_id, start, end, class[, scov]; strata: one label per molecule
    (lengths order) or None; paf: linkstat.loci.read_paf(..., min_block=1000) of the carrier overlaps, or loci: a
    unit-locus function / "positions" (simulation); exclude: molecule ids to drop (chimeras)."""
    ln = lengths.reset_index(drop=True)
    st = None if strata is None else np.asarray(strata, dtype=object)
    excluded = 0
    if exclude is not None:
        keep = ~ln.molecule_id.isin(set(exclude))
        excluded = int((~keep).sum())
        if st is not None:
            st = st[keep.to_numpy()]
        ln = ln[keep].reset_index(drop=True)
    features = features[features.molecule_id.isin(ln.molecule_id)]
    if loci is None and paf is not None:
        loci = loci_from_paf(paf)
    st_by_mol = None if st is None else pd.Series(st, index=ln.molecule_id.to_numpy())
    rows = []
    for a, b in pairs:
        fa, fb = features[features["class"] == a], features[features["class"] == b]
        row = dict(class_a=a, class_b=b, molecules=len(ln), excluded_molecules=excluded)
        if fa.empty or fb.empty:
            rows.append(dict(row, flags="class empty"))
            continue
        r = net_linked_test(
            ln, fa, fb, strata=st, loci=loci, w=w, units_d=units_d, seed=seed, n_mc=n_mc, same_orf=same_orf
        )
        flags = [r.pop("flags")] if r.get("flags") else []
        row.update(r)
        if st_by_mol is not None:
            share = float((st_by_mol.reindex(fa.molecule_id.unique()) == UNRESOLVED).mean())
            row["a_molecules_unresolved"] = share
            if share > 0:
                flags.append(f"{share:.0%} of A molecules in the unresolved host stratum")
        if oe:
            fa2, fb2, _ = drop_same_orf(fa, fb, same_orf) if same_orf is not None else (fa, fb, 0)
            u = cooccurrence_test_units(ln[["molecule_id", "length"]], fa2, fb2, strata=st, d=units_d or 1000)
            row.update(
                oe_fixed_class=a if u.fixed == "A" else b,
                oe_n2u_observed=u.n2.observed,
                oe_n2u_expected=u.n2.expected,
                oe_n2u=u.n2.oe,
                oe_n2u_lo=u.n2.oe_lo,
                oe_n2u_hi=u.n2.oe_hi,
                oe_n2u_p=u.n2.p_two_sided,
                oe_n0u=u.n0.oe,
                oe_n0u_p=u.n0.p_two_sided,
                oe_clustering_index_modelled=u.clustering_index_modelled,
            )
        row["flags"] = "; ".join(flags)
        rows.append(row)
    out = pd.DataFrame(rows)
    for t in ("net_w", "net_mol"):
        if f"{t}_p" in out:
            out[f"{t}_q"] = bh(out[f"{t}_p"].to_numpy(float))
    return out
