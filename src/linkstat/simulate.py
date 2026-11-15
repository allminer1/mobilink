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
Loci: replicon * 1,000,000 + gene index (B genes + 500,000), so the reads covering one locus can be identified.
Options for P2 (all off by default; with the defaults the random-number order and the output are unchanged):
  a_islands       A genes also in islands (cassettes; count ~ Poisson(a_islands), size 1 + Geometric(p_a_island),
                  gaps ~ Exp(island_gap)), placed after B and before the single A genes;
  density_sigma   host heterogeneity: per-genome log-normal multiplier (mean 1) on all gene counts, shared by A and B;
  density_cap     clip that multiplier at +- density_cap x density_sigma on the log scale (after the draw, so the
                  random-number order is unchanged). Without it, 2,000 genomes reach about 30x the base density, and
                  in genomes that crowded, non-overlapping placement packs A genes against B genes (a hard-core
                  contact effect: association under no linkage; P2 step 2.7, probes v4-v5);
  plasmids        replicons per genome (Poisson), lengths uniform in plasmid_length, abundance = host abundance x
                  plasmid_copy, gene density plasmid_density x the chromosome's per base; reads from a plasmid carry
                  the host as `genome` (stratum) and the replicon in `replicon`;
  chimera_rate    (simulate_tables) share of molecules joined from two reads (molecule ids and genome of the first;
                  column `chimera`);
  truth_stratified  also net_strat: anchors per replicon in proportion to its abundance-weighted A genes, anchors
                  touching a B gene redrawn (the definition of the net linked test, amendment 04 draft v2).

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
    a_islands=0.0,
    p_a_island=0.25,
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
            if span >= glen - 1:  # only on short replicons (plasmids); never for 2-6 Mb genomes
                continue
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
    if a_islands > 0:
        n_ai = rng.poisson(a_islands)
        for k in np.maximum(rng.geometric(p_a_island, n_ai), 1) if n_ai else []:
            lens = rng.choice(a_len, k)
            gaps = rng.exponential(island_gap, k).astype(int)
            span = int(lens.sum() + gaps.sum())
            if span >= glen - 1:
                continue
            for _ in range(50):
                s0 = int(rng.integers(0, glen - span))
                e0 = s0 + span
                j = np.searchsorted(occ_s, s0)
                if (j > 0 and occ_e[j - 1] > s0) or (j < len(occ_s) and occ_s[j] < e0):
                    continue
                pos = s0
                for g, gap in zip(lens, gaps, strict=True):
                    j = np.searchsorted(occ_s, pos)
                    occ_s.insert(j, pos)
                    occ_e.insert(j, pos + int(g))
                    A.append((pos, pos + int(g)))
                    pos += int(g) + int(gap)
                break
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
    return_genomes=False,
    a_islands=0.0,
    p_a_island=0.25,
    density_sigma=0.0,
    plasmids=0.0,
    plasmid_length=(20_000, 150_000),
    plasmid_copy=5.0,
    plasmid_density=4.0,
    density_cap=None,
    return_replicons=False,
    return_positions=False,
    truth_stratified=False,
):
    """Returns (reads, truth[, loci]); reads = (lengths, rid, start, end, cls) with cls in {"ARG", "MGE"} and
    truth[(mode, w)] = dict(F_A, F_rnd, net). genome_length (min, max) and abundance_sigma (log-normal sigma) default
    to the step 1.2 community; (L, L) and 0 give equal genomes (the random-number order is the same either way)."""
    rng = np.random.default_rng(seed)
    a_median = float(np.median(a_len))
    glens = rng.integers(int(genome_length[0]), int(genome_length[1]) + 1, n_genomes)
    ab = rng.lognormal(0, abundance_sigma, n_genomes)
    mult = rng.lognormal(-(density_sigma**2) / 2, density_sigma, n_genomes) if density_sigma > 0 else np.ones(n_genomes)
    if density_cap is not None and density_sigma > 0:
        c = -(density_sigma**2) / 2
        mult = np.clip(mult, np.exp(c - density_cap * density_sigma), np.exp(c + density_cap * density_sigma))

    def genome(glen, scale, k):
        return make_genome(
            rng,
            int(glen),
            a_len,
            b_len,
            n_a * scale,
            islands * scale,
            p_island,
            island_gap,
            singletons * scale,
            clustered,
            linked,
            link_gap,
            a_median,
            linked_short,
            linked_long,
            a_islands * scale,
            p_a_island,
        )

    genomes = [genome(gl, mult[k], k) for k, gl in enumerate(glens)]
    host = np.arange(n_genomes)
    if plasmids > 0:
        g_mean = (genome_length[0] + genome_length[1]) / 2
        pl, pab, ph = [], [], []
        for k in range(n_genomes):
            for _ in range(rng.poisson(plasmids)):
                plen = int(rng.integers(int(plasmid_length[0]), int(plasmid_length[1]) + 1))
                genomes.append(genome(plen, mult[k] * plasmid_density * plen / g_mean, k))
                pl.append(plen)
                pab.append(ab[k] * plasmid_copy)
                ph.append(k)
        glens, ab, host = np.r_[glens, pl].astype(np.int64), np.r_[ab, pab], np.r_[host, ph].astype(np.int64)
    n_rep = len(genomes)
    truth = {}
    wa = np.concatenate([np.full(len(A), ab[k]) for k, (A, B) in enumerate(genomes)])
    p_g = ab * glens / (ab * glens).sum()
    gi = rng.choice(n_rep, n_truth_anchors, p=p_g)
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
    if truth_stratified:  # own random stream: the default output is unchanged
        rt = np.random.default_rng(seed + 5_555_555)
        n_a_k = np.array([len(A) for A, _ in genomes], float)
        share = ab * n_a_k / max((ab * n_a_k).sum(), 1e-300)
        for mode in ("start", "end"):
            dA = [side_distances(A, B, mode) for A, B in genomes]
            per = []
            for k, (_A, B) in enumerate(genomes):
                n = int(round(n_truth_anchors * share[k]))
                if n == 0:
                    continue
                got = []
                for _ in range(50):
                    g = rt.choice(a_len, 2 * n)
                    s0 = (rt.random(2 * n) * (glens[k] - g)).astype(np.int64)
                    R = np.column_stack([s0, s0 + g])
                    if len(B):
                        i = np.searchsorted(B[:, 0], R[:, 1], "left")
                        touch = np.zeros(len(R), bool)
                        for dd in (-2, -1, 0):
                            kk = np.clip(i + dd, 0, len(B) - 1)
                            touch |= (B[kk, 0] < R[:, 1]) & (B[kk, 1] > R[:, 0])
                        R = R[~touch]
                    got.append(R)
                    if sum(map(len, got)) >= n:
                        break
                R = np.concatenate(got)[:n]
                R = R[np.argsort(R[:, 0])]
                per.append(side_distances(R, B, mode))
            for w in windows:
                fa = truth[(mode, w)]["F_A"]
                hits = np.concatenate([np.r_[r <= w, lft <= w] for r, lft in per])
                fr = float(hits.mean()) if hits.size else np.nan
                truth[(mode, w)].update(F_rnd_strat=fr, net_strat=1 - (1 - fa) / (1 - fr))
        # unit level (A genes within 1 kb with no B between form a unit; anchors drawn from the unit spans)
        units = []
        for A, B in genomes:
            if len(A) == 0:
                units.append(np.zeros((0, 2), np.int64))
                continue
            gap = A[1:, 0] - A[:-1, 1]
            b_between = np.searchsorted(B[:, 0], A[1:, 0]) > np.searchsorted(B[:, 0], A[:-1, 1]) if len(B) else 0
            new_u = np.r_[True, (gap > 1000) | b_between]
            uid = np.cumsum(new_u) - 1
            U = np.column_stack(
                [
                    np.minimum.reduceat(A[:, 0], np.flatnonzero(new_u)),
                    np.maximum.reduceat(A[:, 1], np.flatnonzero(new_u)),
                ]
            )
            units.append(U[np.argsort(U[:, 0])] if uid.size else U)
        spans = np.concatenate([U[:, 1] - U[:, 0] for U in units]) if units else np.array([], int)
        n_u = np.array([len(U) for U in units], float)
        share_u = ab * n_u / max((ab * n_u).sum(), 1e-300)
        wu = np.concatenate([np.full(len(U), ab[k]) for k, U in enumerate(units)])
        for mode in ("start", "end"):
            dU = [side_distances(U, B, mode) for U, (_A, B) in zip(units, genomes, strict=True)]
            per = []
            for k, (_A, B) in enumerate(genomes):
                n = int(round(n_truth_anchors * share_u[k]))
                if n == 0 or spans.size == 0:
                    continue
                got = []
                for _ in range(50):
                    g = rt.choice(spans, 2 * n)
                    g = g[g < glens[k] - 1]
                    s0 = (rt.random(g.size) * (glens[k] - g)).astype(np.int64)
                    R = np.column_stack([s0, s0 + g])
                    if len(B):
                        i = np.searchsorted(B[:, 0], R[:, 1], "left")
                        touch = np.zeros(len(R), bool)
                        for dd in (-2, -1, 0):
                            kk = np.clip(i + dd, 0, len(B) - 1)
                            touch |= (B[kk, 0] < R[:, 1]) & (B[kk, 1] > R[:, 0])
                        R = R[~touch]
                    got.append(R)
                    if sum(map(len, got)) >= n:
                        break
                R = np.concatenate(got)[:n]
                R = R[np.argsort(R[:, 0])]
                per.append(side_distances(R, B, mode))
            for w in windows:
                fu = (np.concatenate([(r <= w) for r, _ in dU]) * wu).sum() + (
                    np.concatenate([(lft <= w) for _, lft in dU]) * wu
                ).sum()
                fu /= 2 * wu.sum()
                hits = np.concatenate([np.r_[r <= w, lft <= w] for r, lft in per]) if per else np.array([])
                fr = float(hits.mean()) if hits.size else np.nan
                truth[(mode, w)].update(F_units=fu, F_rnd_units=fr, net_units=1 - (1 - fu) / (1 - fr))
    rr = rng if read_seed is None else np.random.default_rng(read_seed)
    n_reads = int(total_bases / read_lengths.mean())
    gr = rr.choice(n_rep, n_reads, p=p_g)
    rl = rr.choice(read_lengths, n_reads)
    lengths, rid, st, en, cl, loc, gen, rep, pos = [], [], [], [], [], [], [], [], []
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
        gen.append(np.full(L.size, host[k], dtype=np.int64))
        rep.append(np.full(L.size, k, dtype=np.int64))
        pos.append(rs)
        base += L.size
    reads = (
        np.concatenate(lengths).astype(np.int64),
        np.concatenate(rid).astype(np.int64),
        np.concatenate(st).astype(np.int64),
        np.concatenate(en).astype(np.int64),
        np.concatenate(cl),
    )
    out = (reads, truth)
    if return_loci:
        out += (np.concatenate(loc).astype(np.int64),)
    if return_genomes:
        out += (np.concatenate(gen),)
    if return_replicons:
        out += (np.concatenate(rep),)
    if return_positions:
        out += (np.concatenate(pos).astype(np.int64),)
    return out


@dataclass
class Simulated:
    lengths: pd.DataFrame  # molecule_id, length, genome (source genome of the read)
    features: pd.DataFrame  # molecule_id, start, end, class, scov, locus
    truth: pd.DataFrame  # mode, w, F_A, F_rnd, net


def simulate_tables(
    read_lengths, a_len, b_len, a_class="ARG", b_class="MGE", chimera_rate=0.0, positions=False, **kw
) -> Simulated:
    """simulate() as the two input tables plus a truth table. The tables do not carry the genes' full lengths, so a
    call touching a read end gets scov = 90 (as a partial alignment would; see stats.netlinked.cut_mask) and every
    other call scov = 100. chimera_rate > 0 joins that share of the output molecules from two reads (own random
    stream, seed + 7,777,777), adds `chimera` to the lengths table and keeps the first read's id and genome.
    positions=True adds `replicon` and `pos` (read start on the replicon) for exact loci
    (linkstat.loci.clusters_from_positions); a chimera keeps the first read's position."""
    plasmids = kw.get("plasmids", 0.0) > 0
    out = simulate(
        np.asarray(read_lengths),
        np.asarray(a_len),
        np.asarray(b_len),
        return_loci=True,
        return_genomes=True,
        return_replicons=plasmids or positions,
        return_positions=positions,
        **kw,
    )
    reads, truth, loci, genome = out[:4]
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
    lengths = pd.DataFrame({"molecule_id": ids, "length": ln, "genome": genome})
    if plasmids or positions:
        lengths["replicon"] = out[4]
    if positions:
        lengths["pos"] = out[5]
    if chimera_rate > 0:
        lengths, feat = _join_chimeras(lengths, feat, chimera_rate, kw.get("seed", 1) + 7_777_777)
    return Simulated(lengths, feat, tr)


def _join_chimeras(lengths: pd.DataFrame, feat: pd.DataFrame, rate: float, seed: int):
    """Join random pairs of reads into one molecule so that `rate` of the output molecules are chimeras."""
    rng = np.random.default_rng(seed)
    n = len(lengths)
    c = int(round(rate * n / (1 + rate)))
    pick = rng.choice(n, 2 * c, replace=False)
    first, second = pick[:c], pick[c:]
    L = lengths.length.to_numpy()
    shift = pd.Series(L[first], index=lengths.molecule_id.to_numpy()[second])
    to = pd.Series(lengths.molecule_id.to_numpy()[first], index=lengths.molecule_id.to_numpy()[second])
    m2 = feat.molecule_id.isin(to.index)
    feat = feat.copy()
    feat.loc[m2, "start"] = feat.loc[m2, "start"] + shift.reindex(feat.loc[m2, "molecule_id"]).to_numpy()
    feat.loc[m2, "end"] = feat.loc[m2, "end"] + shift.reindex(feat.loc[m2, "molecule_id"]).to_numpy()
    feat.loc[m2, "molecule_id"] = to.reindex(feat.loc[m2, "molecule_id"]).to_numpy()
    out = lengths.copy()
    out["chimera"] = False
    out.loc[out.index[first], "length"] = L[first] + L[second]
    out.loc[out.index[first], "chimera"] = True
    out = out.drop(index=out.index[second]).reset_index(drop=True)
    feat = feat.sort_values(["molecule_id", "start"], kind="stable").reset_index(drop=True)
    return out, feat
