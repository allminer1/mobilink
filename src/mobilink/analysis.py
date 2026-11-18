"""`mobilink test`: all tests and effect sizes for a set of class pairs, one row per pair.

Per pair (A:B), after the same-ORF rule (A and B calls overlapping by >= 50% of the shorter are removed):
- molecule level:
  - primary: N2-ev (amendment 03, gate G2 failed -> the pre-registered fallback): exact permutation of the modelled
    class within 1 kb bins of its excluded-volume exposure (and strata), conditional on both margins in every cell;
    O, E, O/E with an exact Poisson interval, p values (two-sided, enrichment, depletion), mid-p, TOST (bound 1.25)
    and the Mantel-Haenszel odds ratio with its interval;
  - N0-ev, the smooth (parametric) approximation of N2-ev: O/E and p. Its single-rate presence curve does not fit a
    clustered class, so it is biased toward depletion when it models one (P1 step 1.15b);
  - the registered length-only exposure null N0 (O/E, p);
  - naive share: molecules with A and B / molecules with A;
- distance level: net linked fraction at w (far-end distance, anchored chance term, gene-length weights) with a
  molecule-bootstrap interval; when no A side has an event within w, the exact upper limit of F_A(w) is given;
  optionally the same after random fragmentation to target N50s;
- Benjamini-Hochberg q values across the pairs (primary two-sided p).
Flags: "pseudo-replication not modelled" (always, until locus clustering exists), "E < 5", "sides < 30",
"modelled class clustered" (clustering index of the non-conditioned class > 2: the nulls model it as independent
genes, so O/E can lean toward depletion; P1 step 1.15b).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .features import drop_same_orf
from .nulls.analytic import clustering_index, cooccurrence_test, cooccurrence_test_ev
from .nulls.binned import binned_test, cell_table, cells, length_bins
from .stats import netlinked as nl
from .stats.equivalence import tost_binned
from .stats.multiple import bh


def _net_linked(lengths, fa, fb, w, n_boot, seed, targets, draws):
    fa = fa.assign(scov=fa["scov"] if "scov" in fa else 100.0)[["molecule_id", "start", "end", "scov"]]
    fb = fb.assign(scov=fb["scov"] if "scov" in fb else 100.0)[["molecule_id", "start", "end", "scov"]]
    if targets:
        row = nl.length_robustness(lengths, fa, fb, w=w, targets=targets, draws=draws, n_boot=n_boot, seed=seed)
        data, scov, _ = nl.prepare(lengths, fa, fb, seed=seed)
    else:
        data, scov, info = nl.prepare(lengths, fa, fb, seed=seed)
        A, R, s = nl.analyse(*data, scov, w)
        lo, hi = nl.bootstrap(A, R, w, n_boot, np.random.default_rng(seed + 1)) if n_boot else (np.nan, np.nan)
        row = dict(w=w, **info, **{f"intact_{k}": v for k, v in s.items()}, intact_lo=lo, intact_hi=hi)
    A, _, _ = nl.analyse(*data, scov, w)
    ev = nl.events_within(A["t"], A["e"], w)
    row["intact_events_w"] = ev
    row["intact_F_A_upper_if_no_event"] = nl.zero_event_upper(int((A["t"] >= w).sum())) if ev == 0 else np.nan
    d = A["t"][A["e"] == 1]
    row["median_event_distance_bp"] = float(np.median(d)) if d.size else np.nan
    return row


def analyse_pairs(
    lengths: pd.DataFrame,
    features: pd.DataFrame,
    pairs,
    strata_col: str | None = None,
    w: int = 2000,
    n_boot: int = 200,
    targets=(),
    draws: int = 5,
    seed: int = nl.SEED,
    min_frac: float = 0.8,
    tost_bound: float = 1.25,
) -> pd.DataFrame:
    """lengths: molecule_id, length[, strata_col]; features: molecule_id, start, end, class[, scov]."""
    strata = lengths[strata_col].to_numpy() if strata_col else None
    L, ids = lengths.length.to_numpy(), lengths.molecule_id
    rows = []
    for a, b in pairs:
        fa, fb, n_dual = drop_same_orf(features[features["class"] == a], features[features["class"] == b])
        row = dict(class_a=a, class_b=b, n_dual=n_dual, calls_A=len(fa), calls_B=len(fb))
        if fa.empty or fb.empty:
            rows.append(dict(row, flags="class empty"))
            continue
        ha, hb = ids.isin(fa.molecule_id).to_numpy(), ids.isin(fb.molecule_id).to_numpy()
        r, cond_on, x = cooccurrence_test_ev(lengths, fa, fb, strata=strata, min_frac=min_frac, tost_bound=None)
        reg = cooccurrence_test(
            L,
            ha,
            hb,
            (fa.end - fa.start).to_numpy(float),
            (fb.end - fb.start).to_numpy(float),
            strata=strata,
            min_frac=min_frac,
            tost_bound=None,
        )
        hx, hy = (hb, ha) if cond_on == "A" else (ha, hb)  # (modelled class, conditioned class)
        n2 = binned_test(L, hx, hy, strata=strata, bin_by=x)
        table = cell_table(cells(length_bins(x), strata), hx, hy)
        tost = tost_binned(table, n2.observed, n2.expected, tost_bound) if tost_bound and n2.expected > 0 else {}
        row.update(
            level="molecule",
            null_model="N2-ev",
            conditioned=b if cond_on == "B" else a,
            molecules=n2.n_molecules,
            molecules_A=int(ha.sum()),
            molecules_B=int(hb.sum()),
            observed=n2.observed,
            expected=n2.expected,
            oe_ratio=n2.oe,
            ci_low=n2.oe_lo,
            ci_high=n2.oe_hi,
            p=n2.p_two_sided,
            p_enriched=n2.p_enriched,
            p_depleted=n2.p_depleted,
            p_mid=n2.p_two_sided_mid,
            p_tost=tost.get("p_tost", np.nan),
            mh_or=n2.mh_odds_ratio,
            mh_or_low=n2.mh_or_lo,
            mh_or_high=n2.mh_or_hi,
            naive_share=n2.observed / ha.sum() if ha.sum() else np.nan,
            n0ev_expected=r.expected,
            n0ev_oe=r.oe,
            n0ev_p=r.p_two_sided,
            registered_expected=reg.expected,
            registered_oe=reg.oe,
            registered_p=reg.p_two_sided,
            modelled=a if cond_on == "B" else b,
            clustering_index_modelled=clustering_index(L, ids, fa if cond_on == "B" else fb, min_frac),
        )
        net = _net_linked(lengths[["molecule_id", "length"]], fa, fb, w, n_boot, seed, targets, draws)
        row.update(
            {
                f"net_{k}" if not k.startswith(("intact", "r")) else k: v
                for k, v in net.items()
                if k not in ("n_dual", "A_calls", "B_calls")
            }
        )
        flags = ["pseudo-replication not modelled"]
        if n2.expected < 5:
            flags.append("E < 5")
        if row["clustering_index_modelled"] > 2:
            flags.append("modelled class clustered (O/E may lean toward depletion)")
        if net.get("intact_sides_w", 0) < 30:
            flags.append("sides < 30")
        row["flags"] = "; ".join(flags)
        rows.append(row)
    out = pd.DataFrame(rows)
    if "p" in out:
        out["q"] = bh(out["p"].to_numpy(float))
    return out
