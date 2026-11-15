"""Poisson-binomial distribution: the number of successes among independent Bernoulli trials with unequal p_i.

Under the read-level exposure null, the number of molecules that carry both A and B is Poisson-binomial with
p_i = p_A(L_i) p_B(L_i).

- pmf(p): exact, by the convolution recursion (O(n^2)); for testing and small n.
- sf(k, p) = P(X >= k) and cdf(k, p) = P(X <= k): exact for n <= max_exact, otherwise the saddlepoint approximation
  with Daniels' (1987) second continuity correction, which keeps a small relative error far into the tails.
Trials with p = 0 are dropped and trials with p = 1 shift k, so they do not cost time or accuracy.
"""

from __future__ import annotations

import numpy as np
from scipy import optimize, special, stats

MAX_EXACT = 5_000


def _clean(p):
    p = np.asarray(p, dtype=float).ravel()
    if np.any((p < 0) | (p > 1)) or np.any(np.isnan(p)):
        raise ValueError("probabilities must lie in [0, 1]")
    ones = int(np.sum(p >= 1.0))
    return p[(p > 0) & (p < 1)], ones


def mean_var(p):
    p = np.asarray(p, dtype=float)
    return float(p.sum()), float((p * (1 - p)).sum())


def pmf(p, max_n: int = 20_000) -> np.ndarray:
    """Exact pmf over 0..n (n = len(p)) by the recursion f_j(x) = f_{j-1}(x) (1 - p_j) + f_{j-1}(x - 1) p_j."""
    p = np.asarray(p, dtype=float).ravel()
    if p.size > max_n:
        raise ValueError(f"exact pmf limited to {max_n} trials")
    f = np.zeros(p.size + 1)
    f[0] = 1.0
    for j, q in enumerate(p, start=1):
        f[1 : j + 1] = f[1 : j + 1] * (1 - q) + f[0:j] * q
        f[0] *= 1 - q
    return f


def _sf_exact(k: int, p) -> float:
    if k <= 0:
        return 1.0
    if k > p.size:
        return 0.0
    f = pmf(p, max_n=max(p.size, 1))
    return float(min(1.0, f[k:].sum()))


def _mills(w):
    """R(w) = (1 - Phi(w)) / phi(w), accurate for large w."""
    return np.sqrt(np.pi / 2) * special.erfcx(w / np.sqrt(2))


def _saddle_upper(k: int, p) -> float:
    """P(X >= k), 1 <= k <= n, by the saddlepoint approximation with Daniels' second continuity correction:
    solve K'(t) = k - 1/2; w = sign(t) sqrt(2 (t (k - 1/2) - K(t))), u = 2 sinh(t / 2) sqrt(K''(t));
    P = 1 - Phi(w) + phi(w) (1/u - 1/w)."""
    n = p.size
    x = k - 0.5
    mu, var = mean_var(p)
    lp, lq = np.log(p), np.log1p(-p)

    def kprime(t):
        return float(np.sum(special.expit(t + lp - lq)))

    def kcgf(t):
        return float(np.sum(np.logaddexp(lq, lp + t)))

    if abs(x - mu) < 1e-8 * max(1.0, np.sqrt(var)):
        return float(stats.norm.sf((x - mu) / np.sqrt(var)))
    lo, hi = -1.0, 1.0
    while kprime(lo) > x:
        lo *= 2
    while kprime(hi) < x and hi < 1e4:
        hi *= 2
    t = optimize.brentq(lambda s: kprime(s) - x, lo, hi, xtol=1e-14, rtol=1e-14, maxiter=500)
    q = special.expit(t + lp - lq)
    k2 = float(np.sum(q * (1 - q)))
    w2 = 2 * (t * x - kcgf(t))
    w = np.sign(t) * np.sqrt(max(w2, 0.0))
    u = 2 * np.sinh(t / 2) * np.sqrt(k2)
    if abs(w) < 1e-4 or abs(u) < 1e-12:
        return float(stats.norm.sf((x - mu) / np.sqrt(var)))
    if w > 0:  # upper tail: P = phi(w) (R(w) + 1/u - 1/w), computed in a form that does not underflow early
        val = np.exp(-0.5 * w * w) / np.sqrt(2 * np.pi) * (_mills(w) + 1 / u - 1 / w)
    else:
        val = stats.norm.sf(w) + stats.norm.pdf(w) * (1 / u - 1 / w)
    return float(min(1.0, max(0.0, val))) if n else 0.0


def sf(k: int, p, method: str = "auto", max_exact: int = MAX_EXACT) -> float:
    """P(X >= k)."""
    p, ones = _clean(p)
    k = int(k) - ones
    if k <= 0:
        return 1.0
    if k > p.size:
        return 0.0
    if method == "exact" or (method == "auto" and p.size <= max_exact):
        return _sf_exact(k, p)
    if method not in ("auto", "saddlepoint"):
        raise ValueError(f"unknown method {method!r}")
    if k == p.size:
        return float(np.exp(np.sum(np.log(p))))
    return _saddle_upper(k, p)


def cdf(k: int, p, method: str = "auto", max_exact: int = MAX_EXACT) -> float:
    """P(X <= k) = P(n - X >= n - k), computed as an upper tail of the complementary trials."""
    p, ones = _clean(p)
    k = int(k) - ones
    if k < 0:
        return 0.0
    if k >= p.size:
        return 1.0
    return sf(p.size - k, 1 - p, method=method, max_exact=max_exact)
