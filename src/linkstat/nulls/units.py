"""Units: the read-level nulls for clustered classes (amendment 04 draft, Part 1; P2 step 2.1).

N0-ev and N2-ev model one class as genes placed independently in the room the other class leaves. Real classes come
in clusters (operons, islands, cassettes; clustering index 3-165 in the G1b runs), and a clustered modelled class
biases O/E toward depletion (G2 memo). Here the modelled class is placed as units instead of genes.

Unit: on a molecule, calls of one class whose gaps are <= d form one unit, spanning from the first start to the last
end; g1, g2 are the lengths of its first and last member (by position).

Room for a modelled unit of span s with edge genes g1, g2 (f = min_frac; a call needs >= f of the gene on the molecule):
  molecule without the fixed class   (L + s - f (g1 + g2))+     (= L - (2f - 1) g for a single gene)
  gap between two fixed units        (gap - s)+                  (the unit must fit whole)
  gap between a fixed unit and an end   (gap - f g)+, averaged over the two orientations (g = g1 or g2 outward)
summed over the gaps and averaged over the unit catalogue. Units overlapping the fixed class are not allowed: a
fixed unit is a solid block (its internal gaps are <= d).
Catalogue: units that end >= d from both molecule ends (complete), weighted by 1 / sum_i (L_i - s - 2d)+, the
relative chance that a unit of span s lies complete on a molecule; if no unit is complete, all units with equal
weights (flagged).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .analytic import cooccurrence_test, exposure_c, fit_density
from .binned import binned_test


def make_units(features: pd.DataFrame, d: int = 1000) -> pd.DataFrame:
    """features: molecule_id, start, end of one class. Returns molecule_id, start, end, n_calls, g1, g2."""
    cols = ["molecule_id", "start", "end", "n_calls", "g1", "g2"]
    if len(features) == 0:
        return pd.DataFrame(columns=cols)
    f = features[["molecule_id", "start", "end"]].sort_values(["molecule_id", "start"], kind="stable")
    f = f.reset_index(drop=True)
    first = (f.molecule_id != f.molecule_id.shift()).to_numpy()
    run_end = f.groupby("molecule_id", sort=False).end.cummax().to_numpy()
    prev = np.r_[-np.inf, run_end[:-1]]
    new = first | (f.start.to_numpy() - prev > d)
    uid = np.cumsum(new) - 1
    f["uid"], f["len"] = uid, f.end - f.start
    g = f.groupby("uid", sort=True)
    last_idx = g.end.idxmax()
    out = pd.DataFrame(
        {
            "molecule_id": g.molecule_id.first(),
            "start": g.start.min(),
            "end": g.end.max(),
            "n_calls": g.size(),
            "g1": g.len.first(),
            "g2": f.loc[last_idx.to_numpy(), "len"].to_numpy(),
        }
    )
    return out.reset_index(drop=True)[cols]


@dataclass
class Catalogue:
    s: np.ndarray  # span
    g1: np.ndarray
    g2: np.ndarray
    w: np.ndarray  # weights, sum 1
    complete: bool  # False if no unit was complete and all units were used


def _sum_hinge(sorted_L: np.ndarray, csum: np.ndarray, t: np.ndarray) -> np.ndarray:
    """sum_i (L_i - t)+ for each t (sorted_L ascending, csum = cumsum with a leading 0)."""
    k = np.searchsorted(sorted_L, t, "right")
    n = sorted_L.size
    return (csum[-1] - csum[k]) - t * (n - k)


def unit_catalogue(units: pd.DataFrame, lengths: pd.DataFrame, d: int = 1000, margin: int | None = None) -> Catalogue:
    """lengths: molecule_id, length. A unit is complete if it ends >= margin (default d) from both molecule ends."""
    m = d if margin is None else margin
    Lmap = pd.Series(lengths.length.to_numpy(float), index=lengths.molecule_id)
    L = Lmap.reindex(units.molecule_id).to_numpy()
    s = (units.end - units.start).to_numpy(float)
    ok = (units.start.to_numpy() >= m) & (L - units.end.to_numpy() >= m) & (units.start.to_numpy() > 0)
    ok &= units.end.to_numpy() < L
    complete = bool(ok.any())
    if not complete:
        ok = np.ones(len(units), bool)
    sL = np.sort(Lmap.to_numpy(float))
    cs = np.r_[0.0, np.cumsum(sL)]
    room = _sum_hinge(sL, cs, s[ok] + (2 * m if complete else 0))
    w = np.where(room > 0, 1 / np.maximum(room, 1e-300), 0.0) if complete else np.ones(ok.sum())
    if w.sum() == 0:
        w = np.ones(ok.sum())
    return Catalogue(
        s=s[ok],
        g1=units.g1.to_numpy(float)[ok],
        g2=units.g2.to_numpy(float)[ok],
        w=w / w.sum(),
        complete=complete,
    )


def unit_room(lengths, molecule_ids, fixed_units: pd.DataFrame, cat: Catalogue, min_frac: float = 0.8) -> np.ndarray:
    """Room for a modelled unit on each molecule, given where the fixed class's units lie (see the module text)."""
    L = np.asarray(lengths, dtype=float)
    f = min_frac
    x = exposure_c(L, f * (cat.g1 + cat.g2) - cat.s, 1.0, cat.w)
    if len(fixed_units) == 0:
        return x
    ids = pd.Index(molecule_ids)
    u = fixed_units[["molecule_id", "start", "end"]].copy()
    u["m"] = ids.get_indexer(u.molecule_id)
    if (u.m < 0).any():
        raise ValueError("fixed units name molecules missing from the lengths table")
    u = u.sort_values(["m", "start"], kind="stable")
    m, s, e = u.m.to_numpy(), u.start.to_numpy(float), u.end.to_numpy(float)
    first = np.r_[True, m[1:] != m[:-1]]
    last = np.r_[m[1:] != m[:-1], True]
    prev_end = np.r_[0.0, e[:-1]]
    interior = np.clip(s - prev_end, 0, None)[~first]
    head = np.clip(s[first], 0, None)
    tail = np.clip(L[m[last]] - e[last], 0, None)
    with_f = m[first]
    xr = np.zeros(L.size)
    np.add.at(xr, m[~first], exposure_c(interior, cat.s, 1.0, cat.w))

    def end_room(gap):
        return 0.5 * (exposure_c(gap, cat.g1, f, cat.w) + exposure_c(gap, cat.g2, f, cat.w))

    np.add.at(xr, with_f, end_room(head) + end_room(tail))
    x[with_f] = xr[with_f]
    return x


def unit_clustering_index(lengths, molecule_ids, units: pd.DataFrame, cat: Catalogue, min_frac: float = 0.8) -> float:
    """As linkstat.nulls.analytic.clustering_index, counting units instead of calls, with the unit exposure."""
    if len(units) == 0:
        return np.nan
    n = units.groupby("molecule_id").size().reindex(molecule_ids, fill_value=0).to_numpy()
    x = exposure_c(np.asarray(lengths, float), min_frac * (cat.g1 + cat.g2) - cat.s, 1.0, cat.w)
    lam = fit_density(x, int((n > 0).sum()))
    if not np.isfinite(lam) or lam == 0:
        return np.nan
    mu = lam * x
    p1 = -np.expm1(-mu)
    p2 = p1 - mu * np.exp(-mu)
    expected = p2.sum() / p1.sum()
    return float(((n >= 2).sum() / (n > 0).sum()) / expected) if expected > 0 else np.nan


@dataclass
class UnitTest:
    fixed: str  # "A" or "B"
    room: np.ndarray  # room of the modelled class per molecule
    n0: object  # CooccurrenceResult (N0-ev-u)
    n2: object  # BinnedResult (N2-ev-u)
    units_a: int
    units_b: int
    catalogue_complete: bool
    clustering_index_modelled: float


def cooccurrence_test_units(
    lengths: pd.DataFrame,
    fa: pd.DataFrame,
    fb: pd.DataFrame,
    strata=None,
    d: int = 1000,
    min_frac: float = 0.8,
    fixed: str | None = None,
    margin: int | None = None,
    width: int = 1000,
    **kw,
) -> UnitTest:
    """N0-ev-u and N2-ev-u for classes A and B (feature tables molecule_id, start, end).
    fixed: "A", "B", None (the class carried by more molecules, as in amendment 03) or "longer" (the class with the
    longer units). margin: catalogue completeness margin (default d). width: N2-ev-u bin width. kw go to
    cooccurrence_test."""
    ids = lengths.molecule_id
    L = lengths.length.to_numpy()
    has_a, has_b = ids.isin(fa.molecule_id).to_numpy(), ids.isin(fb.molecule_id).to_numpy()
    ua, ub = make_units(fa, d), make_units(fb, d)
    if fixed is None:  # amendment 03 rule: the class carried by more molecules
        fixed = "A" if has_a.sum() > has_b.sum() else "B"
    elif fixed == "longer":  # the class with the longer units (mean observed span); its spans are the hardest to learn
        span_a = float((ua.end - ua.start).mean()) if len(ua) else 0.0
        span_b = float((ub.end - ub.start).mean()) if len(ub) else 0.0
        fixed = "A" if span_a > span_b else "B"
    if fixed == "A":
        cat = unit_catalogue(ub, lengths, d, margin)
        x = unit_room(L, ids, ua, cat, min_frac)
        h_mod, h_fix, u_mod = has_b, has_a, ub
    else:
        cat = unit_catalogue(ua, lengths, d, margin)
        x = unit_room(L, ids, ub, cat, min_frac)
        h_mod, h_fix, u_mod = has_a, has_b, ua
    la = (fa.end - fa.start).to_numpy(float)
    lb = (fb.end - fb.start).to_numpy(float)
    len_mod, len_fix = (lb, la) if fixed == "A" else (la, lb)
    n0 = cooccurrence_test(L, h_mod, h_fix, len_mod, len_fix, strata=strata, min_frac=min_frac, exposure_a=x, **kw)
    n2 = binned_test(L, h_mod, h_fix, strata=strata, bin_by=x, width=width)
    return UnitTest(
        fixed=fixed,
        room=x,
        n0=n0,
        n2=n2,
        units_a=len(ua),
        units_b=len(ub),
        catalogue_complete=cat.complete,
        clustering_index_modelled=unit_clustering_index(L, ids, u_mod, cat, min_frac),
    )
