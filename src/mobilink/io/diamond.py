"""DIAMOND blastx tabular output (-f 6 qseqid sseqid pident length qstart qend sstart send slen qlen evalue bitscore).

Coordinates: DIAMOND qstart/qend are 1-based and qstart > qend on the minus strand. Output coordinates are 0-based,
end exclusive. Port of P0_02c_features.read_hits / collapse (frozen in amendment 02 Part B).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

COLS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qstart",
    "qend",
    "sstart",
    "send",
    "slen",
    "qlen",
    "evalue",
    "bitscore",
]
MIN_SCOV = 80.0


def read_hits(path, min_id: float, min_scov: float = MIN_SCOV, chunksize: int = 5_000_000) -> pd.DataFrame:
    """All hits with identity >= min_id and subject cover >= min_scov (subject cover in %, from sstart/send/slen).
    Adds scov, start, end (0-based, end exclusive) and strand. path None returns an empty table."""
    parts = []
    chunks = (
        []
        if path is None
        else pd.read_csv(
            path, sep="\t", names=COLS, header=None, chunksize=chunksize, dtype={"qseqid": str, "sseqid": str}
        )
    )
    for ch in chunks:
        scov = 100.0 * ((ch.send - ch.sstart).abs() + 1) / ch.slen
        keep = (ch.pident >= min_id) & (scov >= min_scov)
        ch = ch.loc[keep].copy()
        ch["scov"] = scov[keep]
        parts.append(ch)
    d = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=COLS + ["scov"])
    d["start"] = np.minimum(d.qstart, d.qend) - 1
    d["end"] = np.maximum(d.qstart, d.qend)
    d["strand"] = np.where(d.qstart <= d.qend, "+", "-")
    return d


def collapse(d: pd.DataFrame, max_overlap: float = 0.5) -> pd.DataFrame:
    """One call per ORF: greedy by bitscore within each query, dropping a hit that overlaps an accepted hit by more than
    `max_overlap` of the shorter of the two. Ties are broken by e-value, then subject ID, so the result does not depend
    on the input order."""
    if d.empty:
        return d
    d = d.sort_values(
        ["qseqid", "bitscore", "evalue", "sseqid"], ascending=[True, False, True, True], kind="stable"
    ).reset_index(drop=True)
    rid = d.qseqid.to_numpy()
    s = d.start.to_numpy()
    e = d.end.to_numpy()
    keep = np.zeros(len(d), dtype=bool)
    acc_s: list = []
    acc_e: list = []
    cur = None
    for i in range(len(d)):
        if rid[i] != cur:
            cur, acc_s, acc_e = rid[i], [], []
        ok = True
        for a, b in zip(acc_s, acc_e, strict=True):
            ov = min(e[i], b) - max(s[i], a)
            if ov > max_overlap * min(e[i] - s[i], b - a):
                ok = False
                break
        if ok:
            keep[i] = True
            acc_s.append(s[i])
            acc_e.append(e[i])
    return d.loc[keep].reset_index(drop=True)
