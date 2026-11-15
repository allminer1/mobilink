"""Feature calling from database hits, and the same-ORF rule used before any pair test.

Threshold rules (identity %, subject cover 80% for all), as registered (docs/preregistration_H1-H4.md):
  primary  ARG 80  MRG 80  MGE 60
  strict   ARG 80  MRG 80  MGE 80
  relaxed  ARG 60  MRG 60  MGE 60
  mgeperm  ARG 80  MRG 80  MGE 40
Per rule and database, hits passing the rule are collapsed to one call per ORF (linkstat.io.diamond.collapse); the
surviving hit labels the ORF. Port of P0_02c_features and P0_04_gonogo.drop_same_orf (frozen in amendment 02 Part B).
"""

from __future__ import annotations

from collections.abc import Callable, Mapping

import numpy as np
import pandas as pd

from .io.diamond import collapse

RULES: dict[str, dict[str, float]] = {
    "primary": {"ARG": 80, "MRG": 80, "MGE": 60},
    "strict": {"ARG": 80, "MRG": 80, "MGE": 80},
    "relaxed": {"ARG": 60, "MRG": 60, "MGE": 60},
    "mgeperm": {"ARG": 80, "MRG": 80, "MGE": 40},
}
OUT_COLS = [
    "molecule_id",
    "start",
    "end",
    "strand",
    "class",
    "gene",
    "label",
    "evidence",
    "intrinsic",
    "pident",
    "scov",
    "bitscore",
    "qlen",
    "db",
]


def check_hit_lengths(hits: pd.DataFrame, lengths: pd.DataFrame) -> tuple[int, int]:
    """(hits whose query is missing from the lengths table, hits whose qlen differs from the table's length)."""
    ln = pd.Series(lengths.length.to_numpy(), index=lengths.molecule_id.to_numpy())
    ql = ln.reindex(hits.qseqid.to_numpy()).to_numpy().astype(float)
    miss = np.isnan(ql)
    return int(miss.sum()), int(np.sum(ql[~miss] != hits.qlen.to_numpy()[~miss]))


def call_features(
    raw: Mapping[str, pd.DataFrame],
    thresholds: Mapping[str, float],
    labellers: Mapping[str, Callable[[pd.DataFrame], pd.DataFrame]],
    order=("ARG", "MRG", "MGE"),
) -> pd.DataFrame:
    """Feature table for one rule. raw[cls]: hits read at (at most) the rule's identity threshold for that class."""
    out = []
    for cls in order:
        if cls not in raw:
            continue
        d = collapse(raw[cls].loc[raw[cls].pident >= thresholds[cls]].copy())
        if d.empty:
            continue
        out.append(labellers[cls](d).rename(columns={"qseqid": "molecule_id"})[OUT_COLS])
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=OUT_COLS)


def summarise(features: pd.DataFrame, **key) -> list[dict]:
    rows = []
    for c, g in features.groupby("class"):
        rows.append(
            dict(
                key,
                cls=c,
                features=len(g),
                molecules=g.molecule_id.nunique(),
                intrinsic=int(g.intrinsic.astype(bool).sum()),
                mean_len=round(float((g.end - g.start).mean()), 1),
            )
        )
    return rows


def drop_same_orf(fa: pd.DataFrame, fb: pd.DataFrame, min_overlap: float = 0.5, id_col: str = "molecule_id"):
    """Remove A and B features that overlap by >= min_overlap of the shorter (one ORF hit by two databases).
    Returns (fa, fb, number of A features removed)."""
    m = fa.reset_index().merge(fb.reset_index(), on=id_col, suffixes=("_a", "_b"))
    if m.empty:
        return fa, fb, 0
    ov = np.minimum(m.end_a, m.end_b) - np.maximum(m.start_a, m.start_b)
    short = np.minimum(m.end_a - m.start_a, m.end_b - m.start_b)
    hit = m.loc[ov >= min_overlap * short]
    return fa.drop(index=hit.index_a.unique()), fb.drop(index=hit.index_b.unique()), int(hit.index_a.nunique())


def covered_bases(features: pd.DataFrame, molecule_ids, id_col: str = "molecule_id") -> np.ndarray:
    """Bases covered by the union of the given features on each molecule (0 for molecules without features), in the
    order of `molecule_ids`. Used for the excluded-volume correction: features cannot overlap other genes, so a
    molecule carrying B has less room for A."""
    out = pd.Series(0.0, index=pd.Index(molecule_ids))
    if len(features) == 0:
        return out.to_numpy()
    f = features[[id_col, "start", "end"]].sort_values([id_col, "start"], kind="stable")
    mid = f[id_col].to_numpy()
    s = f.start.to_numpy(np.int64)
    e = f.end.to_numpy(np.int64)
    new_mol = np.r_[True, mid[1:] != mid[:-1]]
    # running maximum of end within each molecule, to merge overlapping intervals
    run_end = e.copy()
    for i in range(1, len(e)):
        if not new_mol[i]:
            run_end[i] = max(run_end[i - 1], e[i])
    prev_end = np.r_[0, run_end[:-1]]
    # new bases contributed by each interval beyond the running end of the earlier ones
    contrib = np.where(new_mol, e - s, np.where(s >= prev_end, e - s, np.clip(e - prev_end, 0, None)))
    tot = pd.Series(contrib, index=mid).groupby(level=0).sum()
    out.loc[tot.index] = tot.to_numpy(float)
    return out.to_numpy()
