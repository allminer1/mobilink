"""Labels for hits against the three reference databases of the first use case (ARG, MRG, MGE).

CARD (protein homolog model)  ARG; `intrinsic` flags vanR_in_vanO_cl, rpoB2 and HelR.
BacMet2 EXP                   MRG if the entry names a metal or metalloid (class table column is_mrg), else BIOCIDE.
mobileOG-db beatrix-1.6       MGE for the categories in `mge_categories` (default: the four of decision D3), else
                              MGE_std. DIAMOND cuts subject IDs at the first space, so the category comes from a lookup
                              table keyed on the cut ID.
Mobility core (amendment 02 A1): a set of mobileOG subject IDs (cut as DIAMOND reports them).
Port of P0_02c_features (frozen in amendment 02 Part B).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

INTRINSIC = frozenset({"vanR_in_vanO_cl", "rpoB2", "HelR"})
MGE_D3 = frozenset({"integration/excision", "transfer", "replication/recombination/repair", "phage"})


class DatabaseError(ValueError):
    """A hit cannot be labelled from the class or lookup table."""


def label_card(d: pd.DataFrame) -> pd.DataFrame:
    f = d.sseqid.str.split("|", expand=True)
    d["gene"] = f[3]
    d["label"] = f[2]
    d["evidence"] = "CARD_homolog"
    d["intrinsic"] = d.gene.isin(INTRINSIC)
    d["class"] = "ARG"
    d["db"] = "card_homolog"
    return d


def load_bacmet_classes(path) -> pd.DataFrame:
    c = pd.read_csv(
        path, sep="\t", usecols=["BacMet_ID", "Accession", "Gene_name", "metals", "is_mrg"], keep_default_na=False
    )
    c["is_mrg"] = c.is_mrg.astype(str).str.upper().map({"TRUE": True, "FALSE": False})
    if c.is_mrg.isna().any():
        raise DatabaseError("BacMet class table: is_mrg is not TRUE/FALSE everywhere")
    return c


def bacmet_key(sseqid: pd.Series):
    """BacMet2 EXP reuses 5 IDs (BAC0291, BAC0358, BAC0705-0707) for two proteins each, so the key is ID + accession
    (field 4, or field 3 in headers with fewer than 5 fields)."""
    f = sseqid.str.split("|")
    acc = np.where(f.str.len() >= 5, f.str[3], f.str[2])
    return f.str[0], acc


def label_bacmet(d: pd.DataFrame, classes: pd.DataFrame) -> pd.DataFrame:
    d["label"], d["Accession"] = bacmet_key(d.sseqid)
    m = d.merge(
        classes, left_on=["label", "Accession"], right_on=["BacMet_ID", "Accession"], how="left", validate="many_to_one"
    )
    if m.is_mrg.isna().any() or len(m) != len(d):
        raise DatabaseError(f"{int(m.is_mrg.isna().sum())} BacMet hits have no class entry")
    m["gene"] = m.Gene_name
    m["evidence"] = m.metals.fillna("")
    m["intrinsic"] = False
    m["class"] = np.where(m.is_mrg.astype(bool), "MRG", "BIOCIDE")
    m["db"] = "bacmet2_exp"
    return m


def load_mobileog_lookup(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str)


def label_mobileog(d: pd.DataFrame, lookup: pd.DataFrame, mge_categories=MGE_D3, core=None) -> pd.DataFrame:
    """Class MGE for calls in `mge_categories` (D3: MGE vs MGE_std), or, if `core` (a set of subject IDs, e.g. the
    mobility core) is given, for calls whose subject is in the core (MGE vs MGE_other)."""
    m = d.merge(lookup, on="sseqid", how="left")
    if m.category.isna().any():
        raise DatabaseError(f"{int(m.category.isna().sum())} mobileOG hits are missing from the lookup")
    m["label"] = m.category
    m["intrinsic"] = False
    if core is None:
        m["class"] = np.where(m.label.isin(mge_categories), "MGE", "MGE_std")
    else:
        m["class"] = np.where(m.sseqid.isin(core), "MGE", "MGE_other")
    m["db"] = "mobileog_all"
    return m


def load_mobility_core(path) -> frozenset:
    return frozenset(pd.read_csv(path, sep="\t", usecols=["sseqid"]).sseqid)
