"""P3 step 3.4d: the read-level net linked test on the low-coverage libraries of our first soil study (PRJNA1436833).

Libraries: tables/P3/P3_04d_s1_samples.tsv, role "analysed" (27 of 36: the 8 soils that were sequenced again deep are
analysed from their deep data in step 3.4b, and the control library is left out). Annotation: the run script of the
public runs, unchanged (scripts/P3_04d_s1_annotate.sh -> public/EKR1_*/).
Per library and pooled over the libraries (linkstat.pairs.pool_runs: the library is a stratum, anchors never leave
their library and host, unit loci never join reads of two libraries), with the loader and the test of the frozen
stage-1 script (frozen/amendment05/P3_03_confirm.py: registered classes, Kraken2 species strata, unit loci from the
carrier overlaps, net_test_pairs with the default seed and 9,999 draws).
One addition: the reads of the ONT lambda control strand (DCS) are removed first. A read is removed if it has an
alignment of >= 200 bases to the lambda genome (scripts/P3_04d_s1_lambda.sh; this takes the control reads and the
reads that are chimeras of a control strand and a soil molecule). The Kraken2 taxon cannot be used: the database of the
run script assigns the control reads to Escherichia coli (D78).
Variants: "kraken_species" (primary, as the public soil runs) and "no_strata" (sensitivity). No bins exist for these
libraries, so host rule B of the deep soils (D56) cannot be applied.
Our soils are outside the registration; no test between sites is made (vko8soil rule 8).
usage: python -I P3_04d_s1_test.py OUT_PREFIX [--procs 8]
"""
import argparse
import importlib.util
from multiprocessing import Pool

import numpy as np
import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
PAIRS = [("ARG", "MGE"), ("ARG", "MRG"), ("MRG", "MGE")]
LAMBDA_MIN_ALIGNED = 200


def frozen():
    spec = importlib.util.spec_from_file_location("p3_confirm", f"{L}/frozen/amendment05/P3_03_confirm.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def lambda_reads(lib):
    p = pd.read_csv(f"{L}/public/{lib}/{lib}.lambda.paf.gz", sep="\t", header=None, usecols=[0, 1, 2, 3, 10],
                    names=["q", "ql", "qs", "qe", "al"], dtype={"q": str})
    p = p[p.al >= LAMBDA_MIN_ALIGNED]
    cov = ((p.qe - p.qs) / p.ql).groupby(p.q).sum()
    return set(cov.index), set(cov.index[cov < 0.5])


def load(lib):
    """The frozen loader, then the lambda control reads removed from every table."""
    ln, f, st, paf = frozen().load_run(lib)
    lam, chim = lambda_reads(lib)
    assert lam <= set(ln.molecule_id)
    keep = ~ln.molecule_id.isin(lam).to_numpy()
    info = dict(library=lib, reads=len(ln), bases=int(ln.length.sum()), lambda_reads=int((~keep).sum()),
                lambda_bases=int(ln.length[~keep].sum()), lambda_chimeric_reads=len(chim), lambda_reads_with_calls=int(f.molecule_id.isin(lam).sum()))
    ln, st = ln[keep].reset_index(drop=True), np.asarray(st, dtype=object)[keep]
    f = f[~f.molecule_id.isin(lam)].reset_index(drop=True)
    if len(paf):
        paf = paf[~paf.q.isin(lam) & ~paf.t.isin(lam)].reset_index(drop=True)
    info.update(reads_kept=len(ln), bases_kept=int(ln.length.sum()), n50=frozen().n50(ln.length),
                **{f"calls_{c}": int((f["class"] == c).sum()) for c in ("ARG", "MRG", "MGE")},
                strata=len(set(st)), unresolved_share=float(np.mean(st == "unresolved")))
    return ln, f, st, paf, info


def one(args):
    lib, variant = args
    from linkstat.pairs import net_test_pairs

    ln, f, st, paf, info = load(lib)
    t = net_test_pairs(ln, f, PAIRS, strata=st if variant == "kraken_species" else None, paf=paf)
    t.insert(0, "variant", variant)
    t.insert(0, "library", lib)
    return t, info


def pooled(args):
    libs, variant = args
    from linkstat.pairs import net_test_pairs, pool_runs

    data = []
    for lib in libs:
        ln, f, st, paf, _ = load(lib)
        if variant == "no_strata":
            st = np.array(["all"] * len(ln), dtype=object)
        data.append(dict(run=lib, lengths=ln, features=f, strata=st, paf=paf))
    ln, f, st, paf = pool_runs(data)
    t = net_test_pairs(ln, f, PAIRS, strata=st, paf=paf, oe=False)
    t.insert(0, "variant", variant)
    t.insert(0, "library", "pooled")
    t["libraries"] = len(libs)
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--procs", type=int, default=8)
    a = ap.parse_args()
    s = pd.read_csv(f"{L}/tables/P3/P3_04d_s1_samples.tsv", sep="\t", dtype=str)
    libs = s.library[s.role == "analysed"].tolist()
    assert len(libs) == 27
    variants = ["kraken_species", "no_strata"]
    with Pool(a.procs) as p:
        res = p.map(one, [(x, v) for x in libs for v in variants], chunksize=1)
        pool = p.map(pooled, [(libs, v) for v in variants], chunksize=1)
    info = pd.DataFrame([i for (_, i), (x, v) in zip(res, [(x, v) for x in libs for v in variants]) if v == variants[0]])
    info = s[["library", "soil", "locality"]].merge(info, on="library")
    info.to_csv(f"{a.out}.libraries.tsv", sep="\t", index=False)
    t = pd.concat([r for r, _ in res] + pool, ignore_index=True)
    t = s[["library", "soil", "locality"]].merge(t, on="library", how="right")
    t.to_csv(f"{a.out}.linkstat.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 500)
    print(info.to_string(index=False))
    c = ["library", "variant", "class_a", "class_b", "net_w", "net_w_lo", "net_w_hi", "net_w_p", "net_w_O", "net_w_E",
         "net_w_eff", "net_mol", "net_mol_p", "net_mol_O", "net_mol_E", "mol_ratio", "mol_ratio_lo90", "mol_ratio_hi90",
         "hosts_effective"]
    print(t[t.library == "pooled"][c].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
