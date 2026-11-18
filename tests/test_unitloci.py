import numpy as np
import pandas as pd

from mobilink.loci import PAF_COLS
from mobilink.unitloci import parity_components, unit_loci_from_paf, unit_loci_from_positions


def test_parity_components_and_conflict():
    comp, rel, bad = parity_components(5, np.array([0, 1, 3]), np.array([1, 2, 4]), np.array([1, 1, 0]))
    assert comp[0] == comp[1] == comp[2] and comp[3] == comp[4] and comp[0] != comp[3]
    assert rel[0] ^ rel[1] == 1 and rel[1] ^ rel[2] == 1 and rel[0] ^ rel[2] == 0
    assert not bad.any()
    _, _, bad2 = parity_components(3, np.array([0, 1, 0]), np.array([1, 2, 2]), np.array([1, 1, 1]))
    assert bad2.all()  # 0-1 flips, 1-2 flips, so 0-2 must not flip


def _paf(rows):
    return pd.DataFrame(rows, columns=PAF_COLS)


def test_units_map_through_forward_and_reverse_overlaps():
    # read a: unit at 1000-2000. b holds the same locus at 3000-4000 (forward overlap a[0, 4000) -> b[2000, 6000));
    # c holds it reversed at 4500-5500 (a position x -> 6500 - x); d: a homolog (match density 0.1)
    units = pd.DataFrame(
        {"molecule_id": ["a", "b", "c", "d"], "start": [1000, 3000, 4500, 1000], "end": [2000, 4000, 5500, 2000]}
    )
    paf = _paf(
        [
            ("a", 8000, 0, 4000, "+", "b", 9000, 2000, 6000, 3800, 4000, 60),
            ("a", 8000, 0, 4000, "-", "c", 7000, 2500, 6500, 3800, 4000, 60),
            ("a", 8000, 0, 4000, "+", "d", 5000, 0, 4000, 400, 4000, 0),
        ]
    )
    out = unit_loci_from_paf(units, paf)
    assert out.locus[0] == out.locus[1] == out.locus[2] != out.locus[3]
    assert out.parity[0] == out.parity[1] and out.parity[2] != out.parity[0]
    assert not out.conflict.any()


def test_unit_loci_from_positions():
    lab = unit_loci_from_positions(
        np.array([0, 0, 0, 1]), np.array([100, 400, 5000, 100]), np.array([1100, 1300, 6000, 1100])
    )
    assert lab[0] == lab[1] != lab[2] and lab[3] not in (lab[0], lab[2])
