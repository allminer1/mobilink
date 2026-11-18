"""P3 step 3.5 (V8): CPU time and peak memory of `mobilink prep`, `chimeras` and `test` per public run, and the
release check (the released `mobilink test` against the registered per-run tables).

Inputs: tmp/p35/<run>/{prep,chimeras,test}.time (GNU time -v) and results.tsv, written by P3_06_runtime_one.sh.
Outputs: tables/P3/P3_06_runtime.tsv (one row per run and command), tables/P3/P3_06_runtime.summary.tsv (medians and
ranges per command) and tables/P3/P3_06_release_check.tsv (one row per run).
CPU time is user + system seconds of one process; wall time is recorded but the server was busy, so it is not used.
"""
import gzip
import os
import re

import numpy as np
import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
KEY = ["class_a", "class_b"]


def gnu_time(path):
    t = open(path).read()
    g = lambda pat: float(re.search(pat + r": ([\d.]+)", t).group(1))  # noqa: E731
    w = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)", t).group(1).split(":")
    wall = sum(float(x) * 60 ** i for i, x in enumerate(reversed(w)))
    return dict(cpu_s=g("User time \\(seconds\\)") + g("System time \\(seconds\\)"), wall_s=wall,
                peak_gb=g("Maximum resident set size \\(kbytes\\)") / 1e6,
                exit=int(re.search(r"Exit status: (\d+)", t).group(1)))


def rows_gz(path):
    n = 0
    with gzip.open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 24), b""):
            n += chunk.count(b"\n")
    return n


def main():
    runs = [(r, 1, f"{L}/tables/P3/confirm/per_run/{r}.tsv") for r in pd.read_csv(f"{L}/tables/P3/confirm/p3_final.tsv", sep="\t").run_accession]
    runs += [(r, 2, f"{L}/tables/P3/confirm_ext/per_run/{r}.tsv") for r in pd.read_csv(f"{L}/public/screen/p3x_final.tsv", sep="\t").run_accession]
    rt, chk = [], []
    for run, stage, reg in runs:
        O, P = f"{L}/tmp/p35/{run}", f"{L}/public/{run}/{run}"
        if not os.path.exists(f"{O}/DONE"):
            chk.append(dict(run=run, stage=stage, status="not finished"))
            continue
        ln = pd.read_csv(f"{P}.lengths.tsv.gz", sep="\t", usecols=["length"]).length
        size = dict(reads=len(ln), bases=int(ln.sum()),
                    hit_rows=sum(rows_gz(f"{P}.{d}.all.tsv.gz") for d in ("card", "bacmet", "mobileog")),
                    paf_rows=rows_gz(f"{P}.ava.paf.gz"), paf_gz_mb=os.path.getsize(f"{P}.ava.paf.gz") / 1e6)
        m = re.search(r"(\d+) reads judged, (\d+) flagged", open(f"{O}/chimeras.log").read())
        size.update(chimera_judged=int(m.group(1)), chimera_flagged=int(m.group(2)))
        for cmd in ("prep", "chimeras", "test"):
            rt.append(dict(run=run, stage=stage, command=cmd, **size, **gnu_time(f"{O}/{cmd}.time")))
        a = pd.read_csv(f"{O}/results.tsv", sep="\t").set_index(KEY)
        b = pd.read_csv(reg, sep="\t").set_index(KEY)
        cols = [c for c in a.columns if c in b.columns and c != "flags"]
        num = [c for c in cols if pd.api.types.is_numeric_dtype(a[c]) and pd.api.types.is_numeric_dtype(b[c])]
        x, y = a.loc[b.index, num].to_numpy(float), b[num].to_numpy(float)
        same_nan = np.array_equal(np.isnan(x), np.isnan(y))
        d = np.nanmax(np.abs(x - y)) if np.isfinite(x).any() else 0.0
        worst = num[int(np.nanargmax(np.nanmax(np.abs(x - y), axis=0)))] if d > 0 else ""
        chk.append(dict(run=run, stage=stage, status="done", pairs=len(b), columns_compared=len(num), same_missing=same_nan,
                        max_abs_diff=float(d), column_of_max=worst,
                        flags_equal=bool((a.loc[b.index, "flags"].fillna("") == b["flags"].fillna("")).all())))
    rt, chk = pd.DataFrame(rt), pd.DataFrame(chk)
    rt["cpu_s_per_million_reads"] = rt.cpu_s / rt.reads * 1e6
    rt["cpu_s_per_gbp"] = rt.cpu_s / rt.bases * 1e9
    rt.to_csv(f"{L}/tables/P3/P3_06_runtime.tsv", sep="\t", index=False)
    chk.to_csv(f"{L}/tables/P3/P3_06_release_check.tsv", sep="\t", index=False)
    q = lambda s: f"{s.median():.3g} ({s.min():.3g}-{s.max():.3g})"  # noqa: E731
    sm = rt.groupby("command").agg(runs=("run", "size"), cpu_min=("cpu_s", lambda s: q(s / 60)), peak_gb=("peak_gb", q),
                                   cpu_s_per_million_reads=("cpu_s_per_million_reads", q), cpu_s_per_gbp=("cpu_s_per_gbp", q))
    sm.to_csv(f"{L}/tables/P3/P3_06_runtime.summary.tsv", sep="\t")
    pd.set_option("display.width", 250)
    print(sm.to_string())
    print(rt[["reads", "bases", "hit_rows", "paf_rows"]].drop_duplicates().describe().loc[["min", "50%", "max"]].to_string())
    one = rt.drop_duplicates("run")
    print("chimera flags, share of judged reads:", q(one.chimera_flagged / one.chimera_judged), "; share of all reads:",
          q(one.chimera_flagged / one.reads))
    print(chk.status.value_counts().to_dict())
    if "max_abs_diff" in chk:
        print("release check: max |difference|", chk.max_abs_diff.max(), "; missing values equal in", int(chk.same_missing.sum()),
              "of", int((chk.status == "done").sum()), "; flags equal in", int(chk.flags_equal.sum()))
        print(chk[(chk.max_abs_diff > 1e-9) | ~chk.flags_equal.astype(bool) | ~chk.same_missing.astype(bool)].to_string(index=False))


if __name__ == "__main__":
    main()
