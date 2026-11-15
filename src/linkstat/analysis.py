"""`linkstat test`: all tests and effect sizes for a set of class pairs, one row per pair.

Per pair (A:B), after the same-ORF rule (A and B calls overlapping by >= 50% of the shorter are removed):
- molecule level:
  - primary: N0-ev (exposure null with the excluded-volume correction, amendment 03): O, E, O/E with an exact
    Poisson interval, conditional p values (double saddlepoint; two-sided and enrichment), mid-p, TOST (bound 1.25);
  - the registered length-only exposure null N0 (O/E, p), and N2-ev (exact permutation within 1 kb bins of the
    corrected exposure; O/E, p, Mantel-Haenszel odds ratio);
  - naive share: molecules with A and B / molecules with A;
- distance level: net linked fraction at w (far-end distance, anchored chance term, gene-length weights) with a
  molecule-bootstrap interval; when no A side has an event within w, the exact upper limit of F_A(w) is given;
  optionally the same after random fragmentation to target N50s;
- Benjamini-Hochberg q values across the pairs (primary two-sided p).
Flags: "pseudo-replication not modelled" (always, until locus clustering exists), "E < 5", "sides < 30".
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .features import drop_same_orf
from .nulls.analytic import cooccurrence_test, cooccurrence_test_ev
from .nulls.binned import binned_test
from .stats import netlinked as nl
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
        r, cond_on, x = cooccurrence_test_ev(lengths, fa, fb, strata=strata, min_frac=min_frac, tost_bound=tost_bound)
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
        hx, hy = (hb, ha) if cond_on == "A" else (ha, hb)
        n2 = binned_test(L, hx, hy, strata=strata, bin_by=x)
        row.update(
            level="molecule",
            null_model="N0-ev",
            conditioned=b if cond_on == "B" else a,
            molecules=r.n_molecules,
            molecules_A=int(ha.sum()),
            molecules_B=int(hb.sum()),
            observed=r.observed,
            expected=r.expected,
            oe_ratio=r.oe,
            ci_low=r.oe_lo,
            ci_high=r.oe_hi,
            p=r.p_two_sided,
            p_enriched=r.p_enriched,
            p_depleted=r.p_depleted,
            p_mid=r.p_two_sided_mid,
            p_tost=r.p_tost,
            naive_share=r.observed / ha.sum() if ha.sum() else np.nan,
            registered_expected=reg.expected,
            registered_oe=reg.oe,
            registered_p=reg.p_two_sided,
            n2ev_expected=n2.expected,
            n2ev_oe=n2.oe,
            n2ev_p=n2.p_two_sided,
            n2ev_mh_or=n2.mh_odds_ratio,
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
        if r.expected < 5:
            flags.append("E < 5")
        if net.get("intact_sides_w", 0) < 30:
            flags.append("sides < 30")
        row["flags"] = "; ".join(flags)
        rows.append(row)
    out = pd.DataFrame(rows)
    if "p" in out:
        out["q"] = bh(out["p"].to_numpy(float))
    return out
