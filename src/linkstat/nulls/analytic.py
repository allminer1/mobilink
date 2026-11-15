"""N0: the analytic exposure null for molecule-level co-occurrence (length only).

A molecule of length L can carry a feature of length a at about (L - c a)+ positions. With features called when at
least a fraction f of them lies on the molecule (f = 0.8 for DIAMOND's subject-cover rule), the start position can
range over L - (2f - 1) a bases, so c = 2f - 1 (c = 1 for whole features only). For a class with call lengths a_k:
    x(L)   = sum_k pi_k (L - c a_k)+                      (exposure; pi = weights over the observed call lengths)
    p(L)   = 1 - exp(-lambda x(L))                        (probability that the molecule carries the class)
    lambda : sum_i p(L_i) = number of molecules carrying the class   (moment fit, as in the prototype)
Under independence given length, P(A and B on molecule i) = a_i b_i with a_i = p_A(L_i), b_i = p_B(L_i), and the
expected number of molecules with both is E = sum_i a_i b_i. O/E is reported with an exact Poisson interval on O.

Test. The densities are fitted on the same molecules, which fixes sum_i A_i and sum_i B_i. Treating a_i b_i as known
(a plug-in Poisson-binomial test) then overstates the null variance of O - E by up to a factor 1 / ((1 - a)(1 - b))
and the test is conservative; for common classes it almost never rejects (P1 step 1.11, tables/P1_11_*). The default
test is therefore conditional on both margins: P(O >= observed | sum A, sum B) by the double saddlepoint
approximation (linkstat.stats.conditional; with equal lengths it approximates Fisher's exact test). The reported
null variance is that of O - E(lambda_hat) by the delta method for the moment fits:
    O - E(lambda_hat) ~ sum_i d_i,  d_i = A_i B_i - c_a A_i - c_b B_i,
    c_a = sum_i b_i a'_i / sum_i a'_i,  c_b = sum_i a_i b'_i / sum_i b'_i,  a'_i = d a_i / d lambda_A,
with A_i, B_i independent Bernoulli(a_i), Bernoulli(b_i); with equal molecule lengths this is the hypergeometric
variance n a b (1 - a)(1 - b). test="normal" uses this variance with a continuity correction instead of the
saddlepoint. Like any exact-type test of a discrete count, the conditional test is somewhat conservative when
counts are small; the mid-p values `p_*_mid` (P(O > k) + P(O = k) / 2) are reported for that case. The plug-in
Poisson-binomial p values are reported as `p_*_unconditional`.
mode "mean" uses the mean call length and c = 1, which reproduces the prototype (P0) expectation.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd
from scipy import optimize, stats

from ..stats import conditional as cond
from ..stats import poisson_binomial as pb


def exposure(lengths, feat_len, min_frac: float = 0.8, weights=None) -> np.ndarray:
    """x(L) = sum_k pi_k (L - c a_k)+ with c = 2 min_frac - 1 and pi = weights / sum(weights) (equal by default).
    feat_len may be a scalar (one length)."""
    L = np.asarray(lengths, dtype=float)
    a = np.atleast_1d(np.asarray(feat_len, dtype=float))
    c = 2 * min_frac - 1
    pi = np.full(a.size, 1 / a.size) if weights is None else np.asarray(weights, float) / np.sum(weights)
    ca = c * a
    o = np.argsort(ca)
    ca, pi = ca[o], pi[o]
    cw = np.concatenate(([0.0], np.cumsum(pi)))  # weight of the k smallest c*a
    cwa = np.concatenate(([0.0], np.cumsum(pi * ca)))  # weighted sum of the k smallest c*a
    k = np.searchsorted(ca, L, "left")  # c a_k < L contribute
    return L * cw[k] - cwa[k]


def fit_density(x, n_pos: int) -> float:
    """lambda such that sum_i (1 - exp(-lambda x_i)) = n_pos (inf if every exposed molecule carries the class)."""
    x = np.asarray(x, dtype=float)
    if n_pos <= 0 or x.sum() == 0:
        return 0.0
    if n_pos >= np.count_nonzero(x):
        return np.inf

    def f(lam):
        return np.sum(-np.expm1(-lam * x)) - n_pos

    hi = 1e-6
    while f(hi) < 0 and hi < 1e6:
        hi *= 10
    return float(optimize.brentq(f, 0.0, hi, xtol=1e-300, rtol=1e-14))


def presence(x, lam: float) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return -np.expm1(-lam * x) if np.isfinite(lam) else (x > 0).astype(float)


def poisson_interval(obs: int, alpha: float = 0.05):
    """Exact (Garwood) interval for a Poisson mean given one observed count."""
    lo = 0.0 if obs == 0 else stats.chi2.ppf(alpha / 2, 2 * obs) / 2
    hi = stats.chi2.ppf(1 - alpha / 2, 2 * obs + 2) / 2
    return float(lo), float(hi)


def conditional_variance(xa, xb, lam_a: float, lam_b: float) -> float:
    """Null variance of O - E(lambda_hat) for moment-fitted densities (delta method; see the module docstring)."""
    a, b = presence(xa, lam_a), presence(xb, lam_b)
    da = xa * np.exp(-lam_a * xa) if np.isfinite(lam_a) else np.zeros_like(a)
    db = xb * np.exp(-lam_b * xb) if np.isfinite(lam_b) else np.zeros_like(b)
    ca = np.sum(b * da) / np.sum(da) if np.sum(da) > 0 else 0.0
    cb = np.sum(a * db) / np.sum(db) if np.sum(db) > 0 else 0.0
    ab = a * b
    v = ab * (1 - ab) + ca**2 * a * (1 - a) + cb**2 * b * (1 - b) - 2 * ca * ab * (1 - a) - 2 * cb * ab * (1 - b)
    return float(np.sum(v))


@dataclass
class CooccurrenceResult:
    n_molecules: int
    n_strata: int
    molecules_A: int
    molecules_B: int
    observed: int
    expected: float
    var_null: float  # conditional (delta method, summed over strata)
    var_unconditional: float  # plug-in Poisson-binomial
    z: float
    oe: float
    oe_lo: float
    oe_hi: float
    p_enriched: float  # P(O >= observed | margins): double saddlepoint (default) or conditional normal
    p_depleted: float  # P(O <= observed | margins)
    p_two_sided: float  # 2 min(p_enriched, p_depleted), capped at 1
    p_enriched_mid: float  # mid-p: P(O > observed) + P(O = observed) / 2 (closer to nominal size for small counts)
    p_depleted_mid: float
    p_two_sided_mid: float
    p_tost: float  # equivalence test of 1/bound < O/E < bound (small = "at chance")
    tost_bound: float
    p_enriched_unconditional: float  # plug-in Poisson-binomial (conservative)
    p_depleted_unconditional: float
    lambda_A: float  # NaN when several strata are fitted
    lambda_B: float
    exposure_mode: str
    excluded_volume: bool
    min_frac: float
    test: str

    def as_dict(self):
        return asdict(self)


def cooccurrence_test(
    lengths,
    has_a,
    has_b,
    len_a,
    len_b,
    strata=None,
    mode: str = "mixture",
    min_frac: float = 0.8,
    alpha: float = 0.05,
    method: str = "auto",
    test: str = "saddlepoint",
    tost_bound: float = 1.25,
    occupied_b=None,
) -> CooccurrenceResult:
    """Molecule-level co-occurrence of classes A and B against the exposure null.
    lengths: per molecule; has_a, has_b: boolean per molecule; len_a, len_b: call lengths of each class (arrays).
    strata: optional label per molecule (sample, taxon bin, GC bin, ...). N1: the densities are fitted per stratum,
    E sums over strata and the test conditions on both margins in every stratum.
    test: "saddlepoint" (default; conditional on the margins, linkstat.stats.conditional) or "normal" (delta-method
    conditional variance, normal approximation with continuity correction). method: for the unconditional
    Poisson-binomial p values. tost_bound: equivalence bound r for TOST on O/E (registered 1.25; None skips it).
    occupied_b: bases covered by B calls on each molecule (linkstat.features.covered_bases). If given, the
    excluded-volume correction is applied: genes cannot overlap, so A can only lie in the room B leaves. The null
    then conditions on where B lies: A's exposure on molecule i is x_A(L_i - occupied_i), lambda_A is fitted on it,
    b_i = 1 for molecules with B and 0 otherwise, and the test conditions on the A margin (per stratum).
    Without it (the registered length-only null), O/E is biased low whenever B calls take up a noticeable part of
    the molecules that carry them (P1 step 1.15 feasibility check)."""
    from ..stats import equivalence

    lengths = np.asarray(lengths)
    has_a, has_b = np.asarray(has_a, bool), np.asarray(has_b, bool)
    codes = np.zeros(lengths.size, np.int64) if strata is None else pd.factorize(np.asarray(strata))[0]
    S = int(codes.max()) + 1 if codes.size else 1
    if mode == "mean":
        xa = exposure(lengths, np.mean(len_a), 1.0)
        xb = exposure(lengths, np.mean(len_b), 1.0)
        min_frac = 1.0
    elif mode == "mixture":
        xa = exposure(lengths, len_a, min_frac)
        xb = exposure(lengths, len_b, min_frac)
    else:
        raise ValueError(f"unknown exposure mode {mode!r}")
    a, b = np.zeros(lengths.size), np.zeros(lengths.size)
    lam_a, lam_b = np.zeros(S), np.zeros(S)
    var_c = 0.0
    if occupied_b is not None:
        free = np.clip(lengths - np.asarray(occupied_b, float), 0, None)
        xa = exposure(free, np.mean(len_a), 1.0) if mode == "mean" else exposure(free, len_a, min_frac)
        b = has_b.astype(float)
        lam_b[:] = np.nan
        for s in range(S):
            m = codes == s
            lam_a[s] = fit_density(xa[m], int(has_a[m].sum()))
            a[m] = presence(xa[m], lam_a[s])
            da = xa[m] * np.exp(-lam_a[s] * xa[m]) if np.isfinite(lam_a[s]) else np.zeros(int(m.sum()))
            ca = np.sum(b[m] * da) / np.sum(da) if np.sum(da) > 0 else 0.0
            var_c += float(np.sum(a[m] * (1 - a[m]) * (b[m] - ca) ** 2))
    else:
        for s in range(S):
            m = codes == s
            lam_a[s], lam_b[s] = fit_density(xa[m], int(has_a[m].sum())), fit_density(xb[m], int(has_b[m].sum()))
            a[m], b[m] = presence(xa[m], lam_a[s]), presence(xb[m], lam_b[s])
            var_c += conditional_variance(xa[m], xb[m], lam_a[s], lam_b[s])
    # a call on a molecule the model gives zero exposure (e.g. a short call on a short read in "mean" mode)
    a[has_a & (a == 0)] = 1e-12
    b[has_b & (b == 0)] = 1e-12
    p = a * b
    obs = int((has_a & has_b).sum())
    mu, var_u = pb.mean_var(p)
    sd = np.sqrt(var_c) if var_c > 0 else np.nan
    n_a, n_b = np.bincount(codes, has_a, S).astype(int), np.bincount(codes, has_b, S).astype(int)
    st = None if S == 1 else codes
    if test == "saddlepoint":
        p_hi = cond.tail(a, b, n_a, n_b, obs, upper=True, strata=st)
        p_lo = cond.tail(a, b, n_a, n_b, obs, upper=False, strata=st)
        mid_hi = (p_hi + cond.tail(a, b, n_a, n_b, obs + 1, upper=True, strata=st)) / 2
        mid_lo = (p_lo + (cond.tail(a, b, n_a, n_b, obs - 1, upper=False, strata=st) if obs > 0 else 0.0)) / 2
    elif test == "normal":
        p_hi = float(stats.norm.sf((obs - mu - 0.5) / sd)) if var_c > 0 else float(obs > mu)
        p_lo = float(stats.norm.cdf((obs - mu + 0.5) / sd)) if var_c > 0 else float(obs < mu)
        mid_hi = float(stats.norm.sf((obs - mu) / sd)) if var_c > 0 else float(obs > mu)
        mid_lo = float(stats.norm.cdf((obs - mu) / sd)) if var_c > 0 else float(obs < mu)
    else:
        raise ValueError(f"unknown test {test!r}")
    p_tost = (
        equivalence.tost_analytic(a, b, n_a, n_b, obs, mu, tost_bound, st)["p_tost"]
        if tost_bound and mu > 0
        else np.nan
    )
    lo, hi = poisson_interval(obs, alpha)
    return CooccurrenceResult(
        n_molecules=int(lengths.size),
        n_strata=S,
        molecules_A=int(has_a.sum()),
        molecules_B=int(has_b.sum()),
        observed=obs,
        expected=mu,
        var_null=var_c,
        var_unconditional=var_u,
        z=(obs - mu) / sd if var_c > 0 else np.nan,
        oe=obs / mu if mu > 0 else np.nan,
        oe_lo=lo / mu if mu > 0 else np.nan,
        oe_hi=hi / mu if mu > 0 else np.nan,
        p_enriched=p_hi,
        p_depleted=p_lo,
        p_two_sided=min(1.0, 2 * min(p_hi, p_lo)),
        p_enriched_mid=mid_hi,
        p_depleted_mid=mid_lo,
        p_two_sided_mid=min(1.0, 2 * min(mid_hi, mid_lo)),
        p_tost=p_tost,
        tost_bound=tost_bound if tost_bound else np.nan,
        p_enriched_unconditional=pb.sf(obs, p, method=method),
        p_depleted_unconditional=pb.cdf(obs, p, method=method),
        lambda_A=float(lam_a[0]) if S == 1 else np.nan,
        lambda_B=float(lam_b[0]) if S == 1 else np.nan,
        exposure_mode=mode,
        excluded_volume=occupied_b is not None,
        min_frac=min_frac,
        test=test,
    )
