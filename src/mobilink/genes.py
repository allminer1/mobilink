"""Per-gene context table and the screen for fixed-length carrier reads (laboratory plasmids).

Both are descriptive. The per-gene values have no interval and are not tested: use them to see which genes carry a
class-level result, not to rank genes.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .features import drop_same_orf
from .stats.nettest import net_w_estimate


def nearest_b(lengths: pd.DataFrame, fa: pd.DataFrame, fb: pd.DataFrame, w: int) -> pd.DataFrame:
    """One row per class A call: whether a whole class B call lies within w of it on either side (distance from the
    edge of the A call to the far end of the B call), how many of its two sides the molecule covers out to w, and
    the label and gene of the nearest B call within w."""
    ln = pd.Series(lengths.length.to_numpy(), index=lengths.molecule_id.to_numpy())
    b_by_m = {m: g for m, g in fb.groupby("molecule_id", sort=False)}
    L = ln.reindex(fa.molecule_id.to_numpy()).to_numpy()
    rows = []
    for m, s, e, mol_len in zip(fa.molecule_id.to_numpy(), fa.start.to_numpy(), fa.end.to_numpy(), L, strict=True):
        linked, dist, label, gene = False, np.nan, None, None
        g = b_by_m.get(m)
        if g is not None:
            bs, be = g.start.to_numpy(), g.end.to_numpy()
            d = np.where(bs >= e, be - e, np.where(be <= s, s - bs, np.inf))  # far-end distance; overlaps skipped
            if len(d) and np.isfinite(d).any():
                i = int(np.argmin(d))
                if d[i] <= w:
                    linked, dist = True, float(d[i])
                    label = g.label.iloc[i] if "label" in g else None
                    gene = g.gene.iloc[i] if "gene" in g else None
        rows.append((linked, dist, label, gene, int(s >= w) + int(mol_len - e >= w)))
    out = pd.DataFrame(rows, columns=["linked", "distance", "b_label", "b_gene", "sides_reaching_w"], index=fa.index)
    return out


def gene_table(lengths: pd.DataFrame, features: pd.DataFrame, a: str, b: str, strata=None, w: int = 2000,
               units_d: int = 1000, min_calls: int = 30, seed: int = 20261114, annotation: pd.DataFrame | None = None,
               ) -> pd.DataFrame:
    """One row per gene of class a: calls, molecules, calls with a whole class b call within w, the commonest kind of
    b next to it, and net_w for that gene alone where at least min_calls calls and 30 unit sides reach w.
    annotation: optional table keyed on `label` (e.g. CARD aro_index: mechanism, drug class) joined to the rows."""
    fa, fb = features[features["class"] == a], features[features["class"] == b]
    cols = ["class_a", "class_b", "gene", "label", "calls", "molecules", "calls_with_b", "share_with_b",
            "sides_reaching_w", "net_w", "net_w_units", "top_b_label", "top_b_label_share", "top_b_gene"]
    if fa.empty or fb.empty:
        return pd.DataFrame(columns=cols)
    fa, fb, _ = drop_same_orf(fa, fb, 0.5)
    if fa.empty:
        return pd.DataFrame(columns=cols)
    fa = fa.reset_index(drop=True)
    near = nearest_b(lengths, fa, fb, w)
    fa = pd.concat([fa, near], axis=1)
    rows = []
    for (gene, label), g in fa.groupby(["gene", "label"], sort=False, dropna=False):
        lk = g[g.linked]
        top_label = lk.b_label.value_counts()
        top_gene = lk.b_gene.value_counts()
        net, units = np.nan, np.nan
        if len(g) >= min_calls:
            try:
                r = net_w_estimate(lengths, g, fb, strata=strata, w=w, units_d=units_d, seed=seed)
                net, units = r["net_w"], r["units"]
            except ValueError:
                pass
        rows.append(dict(
            class_a=a, class_b=b, gene=gene, label=label, calls=len(g), molecules=g.molecule_id.nunique(),
            calls_with_b=int(g.linked.sum()), share_with_b=float(g.linked.mean()),
            sides_reaching_w=int(g.sides_reaching_w.sum()), net_w=net, net_w_units=units,
            top_b_label=top_label.index[0] if len(top_label) else None,
            top_b_label_share=float(top_label.iloc[0] / len(lk)) if len(lk) else np.nan,
            top_b_gene=top_gene.index[0] if len(top_gene) else None,
        ))
    t = pd.DataFrame(rows, columns=cols).sort_values(["calls_with_b", "calls"], ascending=False).reset_index(drop=True)
    if annotation is not None and not t.empty:
        t = t.merge(annotation, on="label", how="left")
    return t


def fixed_length_clusters(lengths: pd.DataFrame, features: pd.DataFrame, cls: str = "ARG", min_reads: int = 30,
                          max_median: int = 20_000, tol: float = 0.03, min_share: float = 0.25) -> pd.DataFrame:
    """Genes whose carrier reads pile up at one read length: the signature of a small plasmid sequenced whole, which
    is natural in some samples and a laboratory construct (expression or cloning vector) in others. Rule of the
    survey screen: at least min_reads carrier reads, a median carrier length below max_median, and at least
    min_share of the carriers within tol of that median. One row per gene with enough carriers; `flag` marks the
    clusters. A flagged cluster is a reason to look at the reads (for example against NCBI UniVec), not a verdict."""
    f = features[features["class"] == cls]
    cols = ["gene", "carriers", "median_length", "share_near_median", "flag"]
    if f.empty:
        return pd.DataFrame(columns=cols)
    ln = pd.Series(lengths.length.to_numpy(), index=lengths.molecule_id.to_numpy())
    c = f[["gene", "molecule_id"]].drop_duplicates()
    c["length"] = ln.reindex(c.molecule_id.to_numpy()).to_numpy()
    rows = []
    for gene, g in c.groupby("gene", sort=False):
        if len(g) < min_reads:
            continue
        med = float(g.length.median())
        share = float((np.abs(g.length - med) <= tol * med).mean())
        rows.append((gene, len(g), med, share, bool(med < max_median and share >= min_share)))
    return pd.DataFrame(rows, columns=cols).sort_values("share_near_median", ascending=False).reset_index(drop=True)
