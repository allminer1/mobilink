"""The two core input tables.

lengths:  molecule_id, length                       (one row per read or contig)
features: molecule_id, start, end, class, ...       (0-based, end exclusive; extra annotation columns are kept)

`read_id` is accepted as an alias of `molecule_id`, since read-level tables from earlier steps use it.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

ID = "molecule_id"
ALIASES = {"read_id": ID, "contig_id": ID}


class TableError(ValueError):
    """An input table does not follow the schema."""


def _normalise_id(df: pd.DataFrame, what: str) -> pd.DataFrame:
    if ID not in df.columns:
        for alias in ALIASES:
            if alias in df.columns:
                return df.rename(columns={alias: ID})
        raise TableError(f"{what}: no '{ID}' (or read_id / contig_id) column")
    return df


def read_lengths(path, drop_empty: bool = True) -> pd.DataFrame:
    """Read a lengths table. Molecule IDs must be unique; lengths must be non-negative integers.
    Zero-length molecules are dropped by default (they cannot carry a feature)."""
    d = _normalise_id(pd.read_csv(path, sep="\t", dtype={ID: str, "read_id": str, "contig_id": str}), "lengths")
    if "length" not in d.columns:
        raise TableError("lengths: no 'length' column")
    if not d[ID].is_unique:
        raise TableError("lengths: molecule IDs are not unique")
    if not pd.api.types.is_integer_dtype(d.length) or (d.length < 0).any():
        raise TableError("lengths: 'length' must hold non-negative integers")
    if drop_empty:
        d = d[d.length > 0]
    return d.reset_index(drop=True)


def read_features(path) -> pd.DataFrame:
    d = _normalise_id(pd.read_csv(path, sep="\t", dtype={ID: str, "read_id": str, "contig_id": str}), "features")
    missing = {"start", "end", "class"} - set(d.columns)
    if missing:
        raise TableError(f"features: missing columns {sorted(missing)}")
    return d


def validate_features(features: pd.DataFrame, lengths: pd.DataFrame) -> None:
    """Check coordinates against the lengths table: 0 <= start < end <= length, and every molecule is known."""
    ln = pd.Series(lengths.length.to_numpy(), index=lengths[ID].to_numpy())
    L = ln.reindex(features[ID].to_numpy()).to_numpy()
    if np.isnan(L.astype(float)).any():
        n = int(np.isnan(L.astype(float)).sum())
        raise TableError(f"features: {n} rows name a molecule that is not in the lengths table")
    s, e = features.start.to_numpy(), features.end.to_numpy()
    bad = (s < 0) | (e <= s) | (e > L)
    if bad.any():
        raise TableError(f"features: {int(bad.sum())} rows have coordinates outside 0 <= start < end <= length")


def write_features(features: pd.DataFrame, path) -> None:
    features.to_csv(path, sep="\t", index=False, compression="infer")
