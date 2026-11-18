"""Host strata per molecule from Kraken2 per-read labels (D40).

Rule (rank = "S" by default; "G" as sensitivity):
  - a read classified at the target rank or below gets its ancestor at that rank ("S:562");
  - a read classified above the target rank, but at family level or below, keeps its own taxon as a coarser stratum
    ("G:570", "F:543");
  - unclassified reads and reads classified above family (order, class, phylum, domain, root) form one stratum,
    "unresolved". It mixes genomes, so host co-variation is not removed inside it; the share of carriers in it is
    reported and flagged by the caller.
Taxonomy: Kraken2 `ktaxonomy.tsv` (taxid | parent | rank | depth | name), or names/nodes in that layout.
Kraken2 settings (D43): confidence 0, with at least MIN_HITS k-mers with a database hit per read (else "unresolved").
Kraken2's confidence divides by all k-mers, so at a read accuracy near 90% (about 2.5% of 35-mers error free) no read
reaches 0.1; at confidence 0 the gut run SRR23278721 had 95% of carriers classified, at 0.1 only 0.3%.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

UNRESOLVED = "unresolved"
MIN_HITS = 10
_ORDER = {"R": 0, "D": 1, "K": 2, "P": 3, "C": 4, "O": 5, "F": 6, "G": 7, "S": 8}


def read_taxonomy(path: str) -> pd.DataFrame:
    """Kraken2 ktaxonomy.tsv -> DataFrame indexed by taxid with columns parent, rank."""
    t = pd.read_csv(
        path,
        sep="\t",
        header=None,
        usecols=[0, 2, 4],
        names=["taxid", "parent", "rank"],
        dtype={0: np.int64, 2: np.int64, 4: str},
    )
    t["rank"] = t["rank"].str.strip()
    return t.set_index("taxid")


def ancestor_at(taxonomy: pd.DataFrame, rank: str = "S") -> pd.Series:
    """For every taxid, its ancestor with exactly this rank code (itself if it has it), or 0 if none."""
    ids = taxonomy.index.to_numpy()
    pos = pd.Series(np.arange(ids.size), index=ids)
    parent = pos.reindex(taxonomy.parent.to_numpy()).fillna(-1).to_numpy(np.int64)
    is_rank = (taxonomy["rank"] == rank).to_numpy()
    cur = np.arange(ids.size)
    out = np.zeros(ids.size, np.int64)
    todo = np.ones(ids.size, bool)
    for _ in range(128):
        hit = todo & (cur >= 0) & is_rank[np.maximum(cur, 0)]
        out[hit] = ids[cur[hit]]
        todo &= ~hit
        nxt = np.where(cur >= 0, parent[np.maximum(cur, 0)], -1)
        stop = (nxt == cur) | (nxt < 0)
        todo &= ~stop
        cur = np.where(stop, -1, nxt)
        if not todo.any():
            break
    return pd.Series(out, index=ids)


def host_strata(taxids, taxonomy: pd.DataFrame, rank: str = "S") -> np.ndarray:
    """One stratum label per read from its Kraken2 taxid (0 = unclassified); see the module rule."""
    taxids = np.asarray(taxids, dtype=np.int64)
    anc = ancestor_at(taxonomy, rank).reindex(taxids).fillna(0).to_numpy(np.int64)
    own = taxonomy["rank"].reindex(taxids).fillna("").to_numpy(str)
    base = np.array([_ORDER.get(r[:1], -1) if r else -1 for r in own])
    lab = np.full(taxids.size, UNRESOLVED, dtype=object)
    at = anc > 0
    lab[at] = [f"{rank}:{a}" for a in anc[at]]
    coarse = ~at & (taxids > 0) & (base >= _ORDER["F"])
    lab[coarse] = [f"{r}:{t}" for r, t in zip(own[coarse], taxids[coarse], strict=True)]
    return lab


def strata_for_molecules(
    molecule_ids,
    kraken: pd.DataFrame,
    taxonomy: pd.DataFrame,
    rank: str = "S",
    taxid_col: str = "taxid",
    min_hits: int | None = MIN_HITS,
) -> np.ndarray:
    """Strata in the order of molecule_ids from a per-read table (read_id, taxid[, hit_kmers]) that labels every
    molecule (anchors are placed on all molecules of a stratum, not only on carriers; D43). Reads with fewer than
    min_hits hit k-mers (if the table has hit_kmers and min_hits is not None) and molecules without a row get
    "unresolved"."""
    k = kraken.drop_duplicates("read_id").set_index("read_id")
    tax = (
        k[taxid_col].where(k["hit_kmers"] >= min_hits, 0) if min_hits is not None and "hit_kmers" in k else k[taxid_col]
    )
    return host_strata(tax.reindex(molecule_ids).fillna(0).to_numpy(np.int64), taxonomy, rank)
