import numpy as np
import pandas as pd

from mobilink.hosts import UNRESOLVED, host_strata, read_taxonomy, strata_for_molecules

# root 1 > Bacteria 2 (D) > phylum 10 (P) > order 20 (O) > family 30 (F) > genus 40 (G) > species 50 (S)
# > strain 60 (S1); species 51 in genus 40
ROWS = [
    (1, 1, "R"),
    (2, 1, "D"),
    (10, 2, "P"),
    (20, 10, "O"),
    (30, 20, "F"),
    (40, 30, "G"),
    (50, 40, "S"),
    (60, 50, "S1"),
    (51, 40, "S"),
]


def _tax(tmp_path):
    p = tmp_path / "ktaxonomy.tsv"
    p.write_text("".join(f"{t}\t|\t{par}\t|\t{r}\t|\t0\t|\tname{t}\n" for t, par, r in ROWS))
    return read_taxonomy(str(p))


def test_strata_rule(tmp_path):
    tax = _tax(tmp_path)
    lab = host_strata([60, 50, 51, 40, 30, 20, 2, 0, 999], tax)
    assert list(lab) == ["S:50", "S:50", "S:51", "G:40", "F:30", UNRESOLVED, UNRESOLVED, UNRESOLVED, UNRESOLVED]
    assert list(host_strata([60, 51, 30], tax, rank="G")) == ["G:40", "G:40", "F:30"]


def test_strata_for_molecules_fills_missing_reads(tmp_path):
    tax = _tax(tmp_path)
    kr = pd.DataFrame({"read_id": ["r1", "r2"], "taxid": [60, 0]})
    lab = strata_for_molecules(np.array(["r2", "r1", "r3"]), kr, tax)
    assert list(lab) == [UNRESOLVED, "S:50", UNRESOLVED]


def test_minimum_hit_kmers(tmp_path):
    tax = _tax(tmp_path)
    kr = pd.DataFrame({"read_id": ["r1", "r2"], "taxid": [60, 51], "hit_kmers": [9, 10], "taxid_c01": [0, 51]})
    assert list(strata_for_molecules(np.array(["r1", "r2"]), kr, tax)) == [UNRESOLVED, "S:51"]
    assert list(strata_for_molecules(np.array(["r1", "r2"]), kr, tax, min_hits=None)) == ["S:50", "S:51"]
    assert list(strata_for_molecules(np.array(["r1", "r2"]), kr, tax, taxid_col="taxid_c01")) == [UNRESOLVED, "S:51"]
