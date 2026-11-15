"""Metagenome simulator with known linkage truth, for calibration and benchmarking.

Community: n_genomes linear genomes (2-6 Mb), log-normal abundances. Per genome:
  B genes: islands (count ~ Poisson(islands); size 1 + Geometric(p_island); gaps ~ Exp(island_gap)) plus singletons
           (Poisson(singletons)); clustered=False places the same number of B genes uniformly.
  A genes: Poisson(n_a); a share `linked` (or, with linked_short / linked_long, a share that depends on whether the A
           gene is shorter or longer than the median A length) is placed next to a random B gene, on either side, with
           a gap ~ Exp(link_gap); the rest uniformly. Genes never overlap (rejection sampling).
  Gene lengths are drawn from the given arrays (A and B separately).
Reads: total_bases sampled over genomes in proportion to abundance x length, read lengths drawn from `read_lengths`,
  uniform start. A gene is called if it lies whole on the read, or (min_frac given) if at least min_frac of it does
  (truncated to the read, as a caller's subject-cover rule would).
Truth (genome level, abundance-weighted over A gene sides; no reads involved), for the distance to the start ('start')
  or to the far end ('end') of the nearest B: F_A(w), F_rnd(w) for random anchors (lengths from the A array, anchors
  overlapping a B by >= 50% dropped) and net(w) = 1 - (1 - F_A) / (1 - F_rnd).
Loci: genome * 1,000,000 + gene index (B genes + 500,000), so the reads covering one locus can be identified.

Port of scripts/P1_02_simulate.py (P1 step 1.2); the random-number order is unchanged, so a seed gives the same data.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


def _place(rng, glen, lens, occupied_s, occupied_e, tries=50):
    """Uniform non-overlapping placement of genes with lengths `lens` given occupied intervals (sorted lists)."""
    out = []
    s_arr, e_arr = list(occupied_s), list(occupied_e)
    for g in lens:
        for _ in range(tries):
            s = int(rng.integers(0, glen - g))
            e = s + g
            k = np.searchsorted(s_arr, s)
            if (k > 0 and e_arr[k - 1] > s) or (k < len(s_arr) and s_arr[k] < e):
                continue
            s_arr.insert(k, s)
            e_arr.insert(k, e)
            out.append((s, e))
            break
    return out, s_arr, e_arr


def make_genome(
    rng,
    glen,
    a_len,
    b_len,
    n_a,
    islands,
    p_island,
    island_gap,
    singletons,
    clustered,
    linked,
    link_gap,
    a_median,
    linked_short=None,
    linked_long=None,
):
    bs, be = [], []
    n_isl = rng.poisson(islands)
    sizes = 1 + rng.geometric(p_island, n_isl) - 1 if n_isl else np.array([], int)
    sizes = np.maximum(sizes, 1)
    n_single = rng.poisson(singletons)
    if clustered:
        for k in sizes:
            lens = rng.choice(b_len, k)
            gaps = rng.exponential(island_gap, k).astype(int)
            span = int(lens.sum() + gaps.sum())
            for _ in range(50):
                s0 = int(rng.integers(0, glen - span))
                e0 = s0 + span
                j = np.searchsorted(bs, s0)
                if (j > 0 and be[j - 1] > s0) or (j < len(bs) and bs[j] < e0):
                    continue
                pos = s0
                for g, gap in zip(lens, gaps, strict=True):
                    j = np.searchsorted(bs, pos)
                    bs.insert(j, pos)
                    be.insert(j, pos + int(g))
                    pos += int(g) + int(gap)
                break
        _, bs, be = _place(rng, glen, rng.choice(b_len, n_single), bs, be)
    else:
        _, bs, be = _place(rng, glen, rng.choice(b_len, int(sizes.sum()) + n_single), bs, be)
    B = np.array(list(zip(bs, be, strict=True)), dtype=np.int64).reshape(-1, 2)
    occ_s, occ_e = list(bs), list(be)
    A = []
    for g in rng.choice(a_len, rng.poisson(n_a)):
        p_link = linked if linked_short is None else (linked_short if g < a_median else linked_long)
        done = False
        if len(B) and rng.random() < p_link:
            for _ in range(20):
                k = int(rng.integers(0, len(B)))
                gap = int(rng.exponential(link_gap))
                s = B[k, 1] + gap if rng.random() < 0.5 else B[k, 0] - gap - g
                e = s + g
                if s < 0 or e > glen:
                    continue
                j = np.searchsorted(occ_s, s)
                if (j > 0 and occ_e[j - 1] > s) or (j < len(occ_s) and occ_s[j] < e):
                    continue
                occ_s.insert(j, s)
                occ_e.insert(j, e)
                A.append((s, e))
                done = True
                break
        if not done:
            pl, occ_s, occ_e = _place(rng, glen, [g], occ_s, occ_e)
            A += pl
    A = np.array(A, dtype=np.int64).reshape(-1, 2)
    return A[np.argsort(A[:, 0])], B[np.argsort(B[:, 0])]


def side_distances(A, B, mode):
    """Genome-level distances from each A (rows) to the nearest B on the right and on the left (inf if none)."""
    bs, be = B[:, 0], B[:, 1]
    o = np.argsort(be)
    be_sorted, bs_by_end = be[o], bs[o]
    i = np.searchsorted(bs, A[:, 1], "left")
    j = np.searchsorted(be_sorted, A[:, 0], "right") - 1
    right = np.full(len(A), np.inf)
    left = np.full(len(A), np.inf)
    ok = i < len(bs)
    if mode == "start":
        right[ok] = bs[i[ok]] - A[ok, 1]
    else:
        sufmin = np.minimum.accumulate(be[::-1])[::-1] if len(be) else be
        right[ok] = sufmin[i[ok]] - A[ok, 1]
    ok = j >= 0
    if mode == "start":
        left[ok] = A[ok, 0] - be_sorted[j[ok]]
    else:
        premax = np.maximum.accumulate(bs_by_end) if len(bs_by_end) else bs_by_end
        left[ok] = A[ok, 0] - premax[j[ok]]
    return right, left


def drop_overlapping(R, B):
    """Drop rows of R overlapping any B by >= 50% of the shorter."""
    if len(B) == 0 or len(R) == 0:
        return R
    i = np.searchsorted(B[:, 0], R[:, 0], "right") - 1
    keep = np.ones(len(R), bool)
    for d in (-1, 0, 1):
        k = np.clip(i + d, 0, len(B) - 1)
        ov = np.minimum(R[:, 1], B[k, 1]) - np.maximum(R[:, 0], B[k, 0])
        short = np.minimum(R[:, 1] - R[:, 0], B[k, 1] - B[k, 0])
        keep &= ~(ov >= 0.5 * short)
    return R[keep]


def simulate(
    read_lengths,
    a_len,
    b_len,
    seed=1,
    n_genomes=60,
    total_bases=1.5e9,
    n_a=12,
    islands=20,
    p_island=0.25,
    island_gap=150,
    singletons=30,
    clustered=True,
    linked=0.3,
    link_gap=400,
    linked_short=None,
    linked_long=None,
    windows=(1000, 2000, 5000),
    n_truth_anchors=200_000,
    min_frac=None,
    read_seed=None,
    return_loci=False,
    genome_length=(2_000_000, 6_000_000),
    abundance_sigma=1.0,
):
    """Returns (reads, truth[, loci]); reads = (lengths, rid, start, end, cls) with cls in {"ARG", "MGE"} and
    truth[(mode, w)] = dict(F_A, F_rnd, net). genome_length (min, max) and abundance_sigma (log-normal sigma) default
    to the step 1.2 community; (L, L) and 0 give equal genomes (the random-number order is the same either way)."""
    rng = np.random.default_rng(seed)
    a_median = float(np.median(a_len))
    glens = rng.integers(int(genome_length[0]), int(genome_length[1]) + 1, n_genomes)
    ab = rng.lognormal(0, abundance_sigma, n_genomes)
    genomes = [
        make_genome(
            rng,
            int(gl),
            a_len,
            b_len,
            n_a,
            islands,
            p_island,
            island_gap,
            singletons,
            clustered,
            linked,
            link_gap,
            a_median,
            linked_short,
            linked_long,
        )
        for gl in glens
    ]
    truth = {}
    wa = np.concatenate([np.full(len(A), ab[k]) for k, (A, B) in enumerate(genomes)])
    p_g = ab * glens / (ab * glens).sum()
    gi = rng.choice(n_genomes, n_truth_anchors, p=p_g)
    for mode in ("start", "end"):
        dA = [side_distances(A, B, mode) for A, B in genomes]
        R_all, dR = [], []
        for k, (_A, B) in enumerate(genomes):
            n = int((gi == k).sum())
            g = rng.choice(a_len, n)
            s = (rng.random(n) * (glens[k] - g)).astype(np.int64)
            R = drop_overlapping(np.column_stack([s, s + g]), B)
            R = R[np.argsort(R[:, 0])]
            dR.append(side_distances(R, B, mode))
            R_all.append(np.full(len(R), ab[k]))
        wr = np.concatenate(R_all)
        for w in windows:
            fa = (np.concatenate([(r <= w) for r, _ in dA]) * wa).sum() + (
                np.concatenate([(lft <= w) for _, lft in dA]) * wa
            ).sum()
            fa /= 2 * wa.sum()
            fr = (np.concatenate([(r <= w) for r, _ in dR]) * wr).sum() + (
                np.concatenate([(lft <= w) for _, lft in dR]) * wr
            ).sum()
            fr /= 2 * wr.sum()
            truth[(mode, w)] = dict(F_A=fa, F_rnd=fr, net=1 - (1 - fa) / (1 - fr))
    rr = rng if read_seed is None else np.random.default_rng(read_seed)
    n_reads = int(total_bases / read_lengths.mean())
    gr = rr.choice(n_genomes, n_reads, p=p_g)
    rl = rr.choice(read_lengths, n_reads)
    lengths, rid, st, en, cl, loc = [], [], [], [], [], []
    base = 0
    for k, (A, B) in enumerate(genomes):
        L = np.minimum(rl[gr == k], glens[k] - 1)
        rs = (rr.random(L.size) * (glens[k] - L)).astype(np.int64)
        re = rs + L
        for G, name, off in ((A, "ARG", 0), (B, "MGE", 500_000)):
            if len(G) == 0:
                continue
            ends_sorted = np.maximum.accumulate(G[:, 1])
            i0 = np.searchsorted(G[:, 0], rs, "left")
            i1 = np.searchsorted(ends_sorted, re, "right")
            cnt = np.clip(i1 - i0, 0, None)
            r_idx = np.repeat(np.arange(L.size), cnt)
            g_idx = (
                np.concatenate([np.arange(a, b) for a, b in zip(i0, i1, strict=True) if b > a])
                if cnt.sum()
                else np.array([], int)
            )
            rid.append(base + r_idx)
            st.append(G[g_idx, 0] - rs[r_idx])
            en.append(G[g_idx, 1] - rs[r_idx])
            cl.append(np.full(r_idx.size, name))
            loc.append(k * 1_000_000 + off + g_idx)
            if min_frac is not None:
                glen = G[:, 1] - G[:, 0]
                j = np.clip(i0 - 1, 0, None)
                ok = (i0 >= 1) & (G[j, 1] > rs) & (G[j, 1] <= re) & ((G[j, 1] - rs) >= min_frac * glen[j])
                r_idx = np.flatnonzero(ok)
                rid.append(base + r_idx)
                st.append(np.zeros(r_idx.size, np.int64))
                en.append(G[j[ok], 1] - rs[ok])
                cl.append(np.full(r_idx.size, name))
                loc.append(k * 1_000_000 + off + j[ok])
                j = np.clip(i1, 0, len(G) - 1)
                ok = (
                    (i1 < len(G))
                    & (G[j, 0] >= rs)
                    & (G[j, 0] < re)
                    & (G[j, 1] > re)
                    & ((re - G[j, 0]) >= min_frac * glen[j])
                )
                r_idx = np.flatnonzero(ok)
                rid.append(base + r_idx)
                st.append(G[j[ok], 0] - rs[ok])
                en.append(L[ok])
                cl.append(np.full(r_idx.size, name))
                loc.append(k * 1_000_000 + off + j[ok])
        lengths.append(L)
        base += L.size
    reads = (
        np.concatenate(lengths).astype(np.int64),
        np.concatenate(rid).astype(np.int64),
        np.concatenate(st).astype(np.int64),
        np.concatenate(en).astype(np.int64),
        np.concatenate(cl),
    )
    if return_loci:
        return reads, truth, np.concatenate(loc).astype(np.int64)
    return reads, truth


@dataclass
class Simulated:
    lengths: pd.DataFrame  # molecule_id, length
    features: pd.DataFrame  # molecule_id, start, end, class, scov, locus
    truth: pd.DataFrame  # mode, w, F_A, F_rnd, net


def simulate_tables(read_lengths, a_len, b_len, a_class="ARG", b_class="MGE", **kw) -> Simulated:
    """simulate() as the two input tables plus a truth table. The tables do not carry the genes' full lengths, so a
    call touching a read end gets scov = 90 (as a partial alignment would; see stats.netlinked.cut_mask) and every
    other call scov = 100."""
    reads, truth, loci = simulate(
        np.asarray(read_lengths), np.asarray(a_len), np.asarray(b_len), return_loci=True, **kw
    )
    ln, rid, st, en, cl = reads
    ids = np.array([f"sim{i}" for i in range(ln.size)])
    cut = (st == 0) | (en == ln[rid])
    feat = pd.DataFrame(
        {
            "molecule_id": ids[rid],
            "start": st,
            "end": en,
            "class": np.where(cl == "ARG", a_class, b_class),
            "scov": np.where(cut, 90.0, 100.0),
            "locus": loci,
        }
    )
    tr = pd.DataFrame([dict(mode=m, w=w, **v) for (m, w), v in truth.items()])
    return Simulated(pd.DataFrame({"molecule_id": ids, "length": ln}), feat, tr)
