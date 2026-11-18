import io

import numpy as np
import pandas as pd
import pytest

from mobilink import databases as db
from mobilink.features import RULES, call_features, check_hit_lengths, drop_same_orf
from mobilink.io.diamond import COLS, collapse, read_hits
from mobilink.io.tables import TableError, read_lengths, validate_features


def hits_tsv(rows):
    return io.StringIO("\n".join("\t".join(map(str, r)) for r in rows) + "\n")


def test_read_hits_coordinates_and_filters():
    rows = [
        # qseqid sseqid pident length qstart qend sstart send slen qlen evalue bitscore
        ["r1", "s1", 90, 100, 11, 310, 1, 100, 100, 1000, 1e-30, 200],  # + strand, scov 100
        ["r1", "s2", 90, 100, 610, 311, 1, 100, 100, 1000, 1e-30, 200],  # - strand
        ["r2", "s3", 50, 100, 1, 300, 1, 100, 100, 1000, 1e-30, 200],  # identity too low
        ["r2", "s4", 90, 50, 1, 150, 1, 50, 100, 1000, 1e-30, 100],  # subject cover 50% < 80%
    ]
    d = read_hits(hits_tsv(rows), 80.0)
    assert list(d.sseqid) == ["s1", "s2"]
    assert list(d.start) == [10, 310] and list(d.end) == [310, 610]
    assert list(d.strand) == ["+", "-"]
    assert np.allclose(d.scov, 100.0)


def frame(rows):
    d = pd.DataFrame(rows, columns=["qseqid", "sseqid", "start", "end", "bitscore", "evalue"])
    return d


def test_collapse_keeps_best_and_breaks_ties_by_evalue_then_id():
    d = frame(
        [
            ["r1", "b", 0, 100, 50, 1e-5],
            ["r1", "a", 10, 110, 50, 1e-5],  # same score and e-value as "b": "a" wins on subject ID
            ["r1", "c", 0, 100, 60, 1e-9],  # best score wins overall
            ["r1", "d", 300, 400, 10, 1e-2],  # no overlap: kept
            ["r2", "e", 0, 100, 5, 1e-1],
        ]
    )
    out = collapse(d)
    assert list(out.sseqid) == ["c", "d", "e"]
    # order of input rows does not matter
    assert list(collapse(d.iloc[::-1]).sseqid) == ["c", "d", "e"]
    # without the best hit, the tie between "a" and "b" (same bitscore and e-value) goes to the smaller subject ID
    no_c = d[d.sseqid != "c"]
    assert list(collapse(no_c).sseqid) == ["a", "d", "e"]
    assert list(collapse(no_c.iloc[::-1]).sseqid) == ["a", "d", "e"]


def test_collapse_overlap_threshold_is_strictly_more_than_half_of_shorter():
    d = frame([["r1", "a", 0, 100, 10, 1e-5], ["r1", "b", 50, 250, 5, 1e-3]])  # overlap 50 = exactly half of 100
    assert list(collapse(d).sseqid) == ["a", "b"]
    d = frame([["r1", "a", 0, 100, 10, 1e-5], ["r1", "b", 49, 250, 5, 1e-3]])  # overlap 51 > 50
    assert list(collapse(d).sseqid) == ["a"]


def test_drop_same_orf_removes_both_calls():
    fa = pd.DataFrame({"molecule_id": ["r1", "r1"], "start": [0, 1000], "end": [900, 1900]})
    fb = pd.DataFrame({"molecule_id": ["r1", "r2"], "start": [400, 0], "end": [1300, 900]})  # 500 of 900 overlap
    a, b, n = drop_same_orf(fa, fb)
    assert n == 1 and list(a.start) == [1000] and list(b.molecule_id) == ["r2"]


def test_bacmet_key_uses_accession_field():
    ids = pd.Series(["BAC0291|abc|tr|Q1|Organism", "BAC0001|xyz|P2|Org"])
    bid, acc = db.bacmet_key(ids)
    assert list(bid) == ["BAC0291", "BAC0001"] and list(acc) == ["Q1", "P2"]


def test_label_mobileog_requires_lookup_entries():
    d = pd.DataFrame({"sseqid": ["x", "y"]})
    lk = pd.DataFrame({"sseqid": ["x"], "category": ["phage"]})
    with pytest.raises(db.DatabaseError):
        db.label_mobileog(d.copy(), lk)
    lk = pd.DataFrame({"sseqid": ["x", "y"], "category": ["phage", "stability/transfer/defense"]})
    assert list(db.label_mobileog(d.copy(), lk)["class"]) == ["MGE", "MGE_std"]


def test_call_features_applies_rule_threshold_per_class():
    raw = pd.DataFrame(
        [
            ["r1", "g|ARO:1|lbl|geneA", 85, 300, 1, 300, 1, 100, 100, 1000, 1e-30, 200],
            ["r1", "g|ARO:2|lbl|geneB", 70, 300, 401, 700, 1, 100, 100, 1000, 1e-30, 150],
        ],
        columns=COLS,
    )
    raw["scov"] = 100.0
    raw["start"] = raw.qstart - 1
    raw["end"] = raw.qend
    raw["strand"] = "+"
    f = call_features({"ARG": raw}, RULES["primary"], {"ARG": db.label_card})
    assert list(f.gene) == ["geneA"]
    f = call_features({"ARG": raw}, RULES["relaxed"], {"ARG": db.label_card})
    assert list(f.gene) == ["geneA", "geneB"]


def test_tables_schema_checks(tmp_path):
    p = tmp_path / "l.tsv"
    p.write_text("read_id\tlength\nr1\t100\nr2\t0\n")
    ln = read_lengths(p)
    assert list(ln.columns) == ["molecule_id", "length"] and list(ln.molecule_id) == ["r1"]
    p.write_text("read_id\tlength\nr1\t100\nr1\t50\n")
    with pytest.raises(TableError):
        read_lengths(p)
    ln = pd.DataFrame({"molecule_id": ["r1"], "length": [100]})
    validate_features(pd.DataFrame({"molecule_id": ["r1"], "start": [0], "end": [100]}), ln)
    with pytest.raises(TableError):
        validate_features(pd.DataFrame({"molecule_id": ["r1"], "start": [0], "end": [101]}), ln)
    with pytest.raises(TableError):
        validate_features(pd.DataFrame({"molecule_id": ["r9"], "start": [0], "end": [10]}), ln)
    h = pd.DataFrame({"qseqid": ["r1", "r9"], "qlen": [100, 10]})
    assert check_hit_lengths(h, ln) == (1, 0)
