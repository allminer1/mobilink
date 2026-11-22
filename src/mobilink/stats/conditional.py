"""Conditional null distribution of the co-occurrence count O = sum_i A_i B_i given the margins, optionally per stratum.

Molecule i carries A with probability a_i and B with probability b_i, independently. Conditioning on the number of
molecules with A and with B (in each stratum, if strata are given) removes the dependence on how the densities were
fitted, as Fisher's exact test does for a 2 x 2 table. With equal a_i and equal b_i within each stratum, the
conditional law of O is a sum of independent hypergeometric laws (the exact stratified, CMH-type test).

- tail(...): P(O >= k | margins) or P(O <= k | margins) by Skovgaard's double saddlepoint approximation with the
  second continuity correction (Butler 2007, Saddlepoint Approximations with Applications, section 4.2):
      w = sgn(s) sqrt(2 [(K(0, t0) - t0'z) - (K(s, t) - s x - t'z)]),   x = k - 1/2 (upper tail),
      u = 2 sinh(s / 2) sqrt(|K''(s, t)| / |K''_tt(0, t0)|),
      P(O >= k | z) = 1 - Phi(w) - phi(w) (1/w - 1/u),
  where K is the joint cumulant generating function of (O, margins), (s, t) solves K'(s, t) = (x, z) and t0 solves
  K'_t(0, t0) = z. `log_tilt` multiplies the odds of the AB outcome by exp(log_tilt) (an alternative with a common
  odds-ratio-type association; used for equivalence tests).
- tilt_for_mean(...): the tilt whose conditional (saddlepoint) mean of O equals a target.
- variance(...): conditional variance of O at the null saddlepoint.
- exact_pmf(...): exact conditional pmf (one stratum) by recursion over molecules, for small problems and validation.
Molecules with a_i = b_i = 0 are ignored. A stratum where no molecule, or every eligible molecule, carries A (or B)
fixes those indicators; its coordinate is then dropped from the saddlepoint system.
"""

from __future__ import annotations

import numpy as np
from scipy import special, stats


class _System:
    """Per-molecule log probabilities of the outcomes 00, 10, 01, 11 and the stratum layout."""

    def __init__(self, a, b, n_a, n_b, strata=None, log_tilt=0.0):
        a, b = np.asarray(a, float).copy(), np.asarray(b, float).copy()
        g = np.zeros(a.size, dtype=np.int64) if strata is None else np.asarray(strata)
        if g.dtype.kind not in "iu":
            raise ValueError("strata must be integer codes 0..S-1")
        S = int(g.max()) + 1 if g.size else 1
        n_a = np.atleast_1d(np.asarray(n_a, dtype=np.int64))
        n_b = np.atleast_1d(np.asarray(n_b, dtype=np.int64))
        if n_a.size != S or n_b.size != S:
            raise ValueError("one margin per stratum is needed")
        keep = (a > 0) | (b > 0)
        a, b, g = a[keep], b[keep], g[keep]
        # per stratum, without a pass over all molecules for each stratum (the cost with thousands of host strata)
        elig_a = np.bincount(g, a > 0, S).astype(np.int64)
        elig_b = np.bincount(g, b > 0, S).astype(np.int64)
        bad = np.flatnonzero((n_a > elig_a) | (n_b > elig_b))
        if bad.size:
            raise ValueError(f"stratum {bad[0]}: more molecules carry a class than the model allows")
        none_a, none_b = n_a == 0, n_b == 0
        all_a, all_b = ~none_a & (n_a == elig_a), ~none_b & (n_b == elig_b)
        a[all_a[g] & (a > 0)] = 1.0
        a[none_a[g]] = 0.0
        b[all_b[g] & (b > 0)] = 1.0
        b[none_b[g]] = 0.0
        free_a, free_b = ~(none_a | all_a), ~(none_b | all_b)
        with np.errstate(divide="ignore"):
            la, l1a, lb, l1b = np.log(a), np.log1p(-a), np.log(b), np.log1p(-b)
        self.lp = np.stack([l1a + l1b, la + l1b, l1a + lb, la + lb + log_tilt], axis=1)
        self.g, self.S = g, S
        self.n_a, self.n_b = n_a, n_b
        self.fa, self.fb = np.flatnonzero(free_a), np.flatnonzero(free_b)
        self.free_t = np.concatenate([1 + self.fa, 1 + S + self.fb])
        # strata free in both margins: their positions among the free A and the free B coordinates
        both = free_a & free_b
        self.ab = np.flatnonzero(both)
        self.ab_a = np.cumsum(free_a)[both] - 1
        self.ab_b = np.cumsum(free_b)[both] - 1
        self.z = np.concatenate([n_a, n_b]).astype(float)
        # fixed contributions to O from strata where both indicators are fixed
        self.hi = int(np.minimum(n_a, n_b).sum())
        sizes = np.bincount(g, minlength=S)
        self.lo = int(np.maximum(0, n_a + n_b - sizes).sum())

    def cgf(self, theta):
        S, g = self.S, self.g
        s, ta, tb = theta[0], theta[1 : 1 + S], theta[1 + S :]
        lw = self.lp + np.stack([np.zeros_like(g, float), ta[g], tb[g], s + ta[g] + tb[g]], axis=1)
        lz = special.logsumexp(lw, axis=1)
        q = np.exp(lw - lz[:, None])
        q11 = q[:, 3]
        ma, mb = q[:, 1] + q11, q[:, 2] + q11
        grad = np.concatenate([[q11.sum()], np.bincount(g, ma, S), np.bincount(g, mb, S)])
        # Hessian over the coordinates that are ever used: s (index 0), then the free A and the free B margins.
        # The full (1 + 2S) x (1 + 2S) matrix is never needed, and with thousands of strata it took gigabytes.
        fa, fb = self.fa, self.fb
        na = fa.size
        h = np.zeros((1 + na + fb.size, 1 + na + fb.size))
        h[0, 0] = np.sum(q11 * (1 - q11))
        c_oa, c_ob = np.bincount(g, q11 * (1 - ma), S), np.bincount(g, q11 * (1 - mb), S)
        h[0, 1 : 1 + na] = h[1 : 1 + na, 0] = c_oa[fa]
        h[0, 1 + na :] = h[1 + na :, 0] = c_ob[fb]
        ia, ib = np.arange(1, 1 + na), np.arange(1 + na, 1 + na + fb.size)
        h[ia, ia] = np.bincount(g, ma * (1 - ma), S)[fa]
        h[ib, ib] = np.bincount(g, mb * (1 - mb), S)[fb]
        c_ab = np.bincount(g, q11 - ma * mb, S)[self.ab]
        h[1 + self.ab_a, 1 + na + self.ab_b] = h[1 + na + self.ab_b, 1 + self.ab_a] = c_ab
        return float(lz.sum()), grad, h

    def solve(self, x=None, theta0=None, max_iter=200):
        """Null saddlepoint (x None: s = 0, solve the free t) or joint saddlepoint at O = x."""
        idx = self.free_t if x is None else np.concatenate([[0], self.free_t])
        hs = slice(1, None) if x is None else slice(None)  # the same coordinates in the reduced Hessian of cgf
        tgt = self.z[idx[idx > 0] - 1] if x is None else np.concatenate([[x], self.z[self.free_t - 1]])
        theta = np.zeros(1 + 2 * self.S) if theta0 is None else np.array(theta0, float)
        k, gr, h = self.cgf(theta)
        f = k - theta[idx] @ tgt
        scale = max(1.0, float(np.max(np.abs(tgt)))) if tgt.size else 1.0
        for _ in range(max_iter):
            r = gr[idx] - tgt
            if r.size == 0 or np.max(np.abs(r)) < 1e-9 * scale:
                return theta, True
            try:
                step = np.linalg.solve(h[hs, hs], r)
            except np.linalg.LinAlgError:
                return theta, False
            lam = 1.0
            while True:
                new = theta.copy()
                new[idx] -= lam * step
                k2, g2, h2 = self.cgf(new)
                f2 = k2 - new[idx] @ tgt
                if f2 <= f + 1e-12 * max(1.0, abs(f)) or lam < 1e-12:
                    break
                lam /= 2
            if abs(f - f2) < 1e-15 * max(1.0, abs(f)) and lam < 1e-12:
                theta, k, gr, h = new, k2, g2, h2
                break
            theta, k, gr, h, f = new, k2, g2, h2, f2
        r = gr[idx] - tgt
        return theta, bool(r.size == 0 or np.max(np.abs(r)) < 1e-6 * scale)


def _schur_var(h):
    """K_ss - K_st K_tt^-1 K_ts over the free coordinates (h: the reduced Hessian of _System.cgf)."""
    if h.shape[0] == 1:
        return float(h[0, 0])
    hst = h[0, 1:]
    return float(h[0, 0] - hst @ np.linalg.solve(h[1:, 1:], hst))


def variance(a, b, n_a, n_b, strata=None) -> float:
    """Conditional variance of O at the null saddlepoint."""
    sy = _System(a, b, n_a, n_b, strata)
    t0, _ = sy.solve()
    return _schur_var(sy.cgf(t0)[2])


def tilt_for_mean(a, b, n_a, n_b, target_mean: float, strata=None) -> float:
    """Tilt s whose joint saddlepoint has O = target_mean given the margins (the conditional mean, to first order)."""
    sy = _System(a, b, n_a, n_b, strata)
    t0, _ = sy.solve()
    th, ok = sy.solve(x=float(target_mean), theta0=t0)
    if not ok:
        raise ValueError("no saddlepoint for this target mean (outside the support?)")
    return float(th[0])


def tail(a, b, n_a, n_b, k: int, upper: bool = True, strata=None, log_tilt: float = 0.0) -> float:
    """P(O >= k | margins) if upper, else P(O <= k | margins), by the double saddlepoint approximation."""
    sy = _System(a, b, n_a, n_b, strata, log_tilt)
    lo_s, hi_s = sy.lo, sy.hi
    if upper:
        if k <= lo_s:
            return 1.0
        if k > hi_s:
            return 0.0
        x = k - 0.5
    else:
        if k >= hi_s:
            return 1.0
        if k < lo_s:
            return 0.0
        x = k + 0.5
    t0, ok0 = sy.solve()
    k0, g0, h0 = sy.cgf(t0)
    cmean, cvar = g0[0], _schur_var(h0)
    th, ok = sy.solve(x=x, theta0=t0)
    s = th[0]
    if not (ok0 and ok) or abs(s) < 1e-6 or cvar <= 0:
        zc = (x - cmean) / np.sqrt(cvar) if cvar > 0 else (np.inf if x < cmean else -np.inf)
        return float(stats.norm.sf(zc) if upper else stats.norm.cdf(zc))
    kh, _, hh = sy.cgf(th)
    zt = sy.z[sy.free_t - 1]
    w2 = 2 * ((k0 - t0[sy.free_t] @ zt) - (kh - s * x - th[sy.free_t] @ zt))
    w = np.sign(s) * np.sqrt(max(w2, 0.0))
    _, ld_joint = np.linalg.slogdet(hh)
    ld_null = np.linalg.slogdet(h0[1:, 1:])[1] if sy.free_t.size else 0.0
    u = 2 * np.sinh(s / 2) * np.exp(0.5 * (ld_joint - ld_null))
    if abs(w) < 1e-8:
        zc = (x - cmean) / np.sqrt(cvar)
        return float(stats.norm.sf(zc) if upper else stats.norm.cdf(zc))
    if upper:
        val = stats.norm.sf(w) - stats.norm.pdf(w) * (1 / w - 1 / u)
    else:
        val = stats.norm.cdf(w) + stats.norm.pdf(w) * (1 / w - 1 / u)
    return float(min(1.0, max(0.0, val)))


def exact_pmf(a, b, n_a: int, n_b: int) -> np.ndarray:
    """Exact P(O = o | S_A = n_a, S_B = n_b), o = 0..min(n_a, n_b), by recursion over molecules
    (array of size (n_a + 1) (n_b + 1) (min + 1); one stratum; for small problems)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = min(n_a, n_b)
    f = np.zeros((n_a + 1, n_b + 1, m + 1))
    f[0, 0, 0] = 1.0
    for ai, bi in zip(a, b, strict=True):
        g = f * ((1 - ai) * (1 - bi))
        g[1:, :, :] += f[:-1, :, :] * (ai * (1 - bi))
        g[:, 1:, :] += f[:, :-1, :] * ((1 - ai) * bi)
        g[1:, 1:, 1:] += f[:-1, :-1, :-1] * (ai * bi)
        f = g
    col = f[n_a, n_b, :]
    return col / col.sum()
