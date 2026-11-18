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
- Chimeras: molecules listed in `exclude` (e.g. flagged by `linkstat chimeras`) are removed before anything else; the
  count is reported.
- Reporting (v0.2): a test with no observed event and fewer than one expected gives no estimate (flagged) instead of
  0 with a zero-width interval; upper interval limits are clipped at 1 (the measure's maximum); fewer than 3
  effective hosts are flagged.
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
    unit-locus function / "positions" (simulation) / "molecule" (contigs: all A units of a molecule are one clump);
    exclude: molecule ids to drop (chimeras)."""
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
        try:
            r = net_linked_test(
                ln, fa, fb, strata=st, loci=loci, w=w, units_d=units_d, seed=seed, n_mc=n_mc, same_orf=same_orf
            )
        except ValueError as e:  # too few molecules or calls for anchors and sides
            rows.append(dict(row, flags=f"test not possible ({e})"))
            continue
        flags = [r.pop("flags")] if r.get("flags") else []
        for t in ("net_w", "net_mol"):
            # reporting (v0.2, P2 step 2.6): no observed event and < 1 expected -> no estimate (it would read 0 +- 0)
            if r.get(f"{t}_O") == 0 and (r.get(f"{t}_E") or 0) < 1:
                r[t] = r[f"{t}_lo"] = r[f"{t}_hi"] = np.nan
                flags.append(f"{t}: no events (E = {r.get(f'{t}_E') or 0:.2g}; no estimate)")
            elif np.isfinite(r.get(f"{t}_hi", np.nan)) and r[f"{t}_hi"] > 1:  # the measure cannot exceed 1
                r[f"{t}_hi"] = 1.0
        if st is not None and r.get("hosts_effective", np.inf) < 3:
            flags.append(f"few effective hosts (G_eff {r['hosts_effective']:.1f}; intervals wide)")
        row.update(r)
        if st_by_mol is not None:
            lab = st_by_mol.reindex(fa.molecule_id.unique()).astype(str)
            share = float(((lab == UNRESOLVED) | lab.str.endswith("|" + UNRESOLVED)).mean())
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


def pool_runs(runs) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, pd.DataFrame | None]:
    """Several runs as one dataset for net_test_pairs (amendment 05 Part C). runs: iterable of dicts with keys run
    (label), lengths (molecule_id, length), features (molecule_id, start, end, class[, scov]), strata (one label per
    molecule, lengths order) and optionally paf (read_paf output). Molecule ids get the prefix "<run>|" and host strata
    become "<run>|<stratum>", so anchors never leave their run and host, and unit loci never join reads of two runs."""
    lns, fs, sts, pafs = [], [], [], []
    for r in runs:
        tag = f"{r['run']}|"
        ln = r["lengths"][["molecule_id", "length"]].copy()
        ln["molecule_id"] = tag + ln.molecule_id.astype(str)
        f = r["features"].copy()
        f["molecule_id"] = tag + f.molecule_id.astype(str)
        st = np.asarray(r["strata"], dtype=object)
        if st.size != len(ln):
            raise ValueError(f"run {r['run']}: {st.size} strata for {len(ln)} molecules")
        lns.append(ln)
        fs.append(f)
        sts.append(np.array([tag + str(x) for x in st], dtype=object))
        if r.get("paf") is not None and len(r["paf"]):
            p = r["paf"].copy()
            p["q"], p["t"] = tag + p.q.astype(str), tag + p.t.astype(str)
            pafs.append(p)
    paf = pd.concat(pafs, ignore_index=True) if pafs else None
    return pd.concat(lns, ignore_index=True), pd.concat(fs, ignore_index=True), np.concatenate(sts), paf

