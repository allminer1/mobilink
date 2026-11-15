"""Conditional null distribution of the co-occurrence count O = sum_i A_i B_i given both margins.

Molecule i carries A with probability a_i and B with probability b_i, independently. Conditioning on
S_A = sum A_i = n_A and S_B = sum B_i = n_B removes the dependence on how the densities were fitted (as Fisher's
exact test does for a 2 x 2 table); with equal a_i and equal b_i the conditional law of O is hypergeometric.

- tail(...): P(O >= k | S_A, S_B) or P(O <= k | S_A, S_B) by Skovgaard's double saddlepoint approximation with the
  second continuity correction (Butler 2007, Saddlepoint Approximations with Applications, section 4.2):
      w = sgn(s) sqrt(2 [(K(0, t0) - t0'z) - (K(s, t) - s x - t'z)]),   x = k - 1/2,
      u = 2 sinh(s / 2) sqrt(|K''(s, t)| / |K''_tt(0, t0)|),
      P(O >= k | z) = 1 - Phi(w) - phi(w) (1/w - 1/u),
  where K is the joint cumulant generating function of (O, S_A, S_B), (s, t) solves K'(s, t) = (x, z) and t0 solves
  K'_t(0, t0) = z.
- exact_pmf(...): exact conditional pmf by recursion over molecules (small problems; used to validate tail()).
- variance(...): conditional variance at the null saddlepoint (Schur complement of K'').
"""

from __future__ import annotations

import numpy as np
from scipy import special, stats

# outcome vectors (O, S_A, S_B) for (A, B) = 00, 10, 01, 11
_V = np.array([[0, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], dtype=float)


def _logp(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    keep = (a > 0) | (b > 0)
    a, b = a[keep], b[keep]
    with np.errstate(divide="ignore"):
        la, l1a, lb, l1b = np.log(a), np.log1p(-a), np.log(b), np.log1p(-b)
    return np.stack([l1a + l1b, la + l1b, l1a + lb, la + lb], axis=1)  # (n, 4)


def _cgf(lp, theta):
    """K(theta), gradient and Hessian for theta = (s, t1, t2)."""
    lw = lp + _V @ theta
    lz = special.logsumexp(lw, axis=1)
    q = np.exp(lw - lz[:, None])
    m = q @ _V  # per-molecule tilted means (n, 3)
    grad = m.sum(0)
    second = np.einsum("nk,ki,kj->ij", q, _V, _V)
    hess = second - m.T @ m
    return float(lz.sum()), grad, hess


def _solve(lp, target, free=(0, 1, 2), theta0=None, max_iter=200):
    """Minimise K(theta) - theta'target over the free coordinates (others fixed at theta0) by damped Newton."""
    theta = np.zeros(3) if theta0 is None else np.array(theta0, float)
    idx = np.array(free)
    tgt = np.asarray(target, float)
    k, g, h = _cgf(lp, theta)
    f = k - theta[idx] @ tgt
    for _ in range(max_iter):
        r = g[idx] - tgt
        hh = h[np.ix_(idx, idx)]
        try:
            step = np.linalg.solve(hh, r)
        except np.linalg.LinAlgError:
            return theta, False
        lam = 1.0
        while True:
            new = theta.copy()
            new[idx] -= lam * step
            k2, g2, h2 = _cgf(lp, new)
            f2 = k2 - new[idx] @ tgt
            if f2 <= f + 1e-12 * max(1.0, abs(f)) or lam < 1e-10:
                break
            lam /= 2
        theta, k, g, h, f_old, f = new, k2, g2, h2, f, f2
        if np.max(np.abs(g[idx] - tgt)) < 1e-9 * max(1.0, np.max(np.abs(tgt))) or abs(f_old - f) < 1e-14:
            return theta, bool(np.max(np.abs(g[idx] - tgt)) < 1e-6 * max(1.0, np.max(np.abs(tgt))))
    return theta, False


def variance(a, b, n_a: int, n_b: int) -> float:
    """Conditional variance of O at the null saddlepoint: K_ss - K_st K_tt^-1 K_ts at (0, t0)."""
    lp = _logp(a, b)
    t0, _ = _solve(lp, [n_a, n_b], free=(1, 2))
    _, _, h = _cgf(lp, t0)
    return float(h[0, 0] - h[0, 1:] @ np.linalg.solve(h[1:, 1:], h[1:, 0]))


def tail(a, b, n_a: int, n_b: int, k: int, upper: bool = True) -> float:
    """P(O >= k | S_A = n_a, S_B = n_b) if upper, else P(O <= k | ...), by the double saddlepoint approximation."""
    lo_support, hi_support = 0, min(n_a, n_b)
    if upper:
        if k <= lo_support:
            return 1.0
        if k > hi_support:
            return 0.0
        x = k - 0.5
    else:
        if k >= hi_support:
            return 1.0
        if k < lo_support:
            return 0.0
        x = k + 0.5
    lp = _logp(a, b)
    z = np.array([n_a, n_b], float)
    t0, ok0 = _solve(lp, z, free=(1, 2))
    k0, _, h0 = _cgf(lp, t0)
    th, ok = _solve(lp, np.array([x, n_a, n_b]), free=(0, 1, 2), theta0=t0)
    kh, _, hh = _cgf(lp, th)
    cvar = h0[0, 0] - h0[0, 1:] @ np.linalg.solve(h0[1:, 1:], h0[1:, 0])
    cmean = _cgf(lp, t0)[1][0]
    s = th[0]
    if not (ok0 and ok) or abs(s) < 1e-6:
        zc = (x - cmean) / np.sqrt(cvar)  # normal approximation with continuity correction
        return float(stats.norm.sf(zc) if upper else stats.norm.cdf(zc))
    w2 = 2 * ((k0 - t0[1:] @ z) - (kh - s * x - th[1:] @ z))
    w = np.sign(s) * np.sqrt(max(w2, 0.0))
    u = 2 * np.sinh(s / 2) * np.sqrt(np.linalg.det(hh) / np.linalg.det(h0[1:, 1:]))
    if abs(w) < 1e-8:
        zc = (x - cmean) / np.sqrt(cvar)
        return float(stats.norm.sf(zc) if upper else stats.norm.cdf(zc))
    p_ge = stats.norm.sf(w) - stats.norm.pdf(w) * (1 / w - 1 / u)  # P(O >= x + 1/2 | z)
    val = p_ge if upper else 1 - p_ge
    if not upper and w < 0:  # compute the lower tail directly so it does not lose precision
        val = stats.norm.cdf(w) + stats.norm.pdf(w) * (1 / w - 1 / u)
    return float(min(1.0, max(0.0, val)))


def exact_pmf(a, b, n_a: int, n_b: int) -> np.ndarray:
    """Exact P(O = o | S_A = n_a, S_B = n_b), o = 0..min(n_a, n_b), by recursion over molecules
    (array of size (n_a + 1) (n_b + 1) (min + 1); for small problems)."""
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
