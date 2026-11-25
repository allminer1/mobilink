"""Power and design: what a given sequencing depth can detect.

Read-level co-occurrence test (N0, conditional):
  For n molecules with lengths drawn from a length distribution, E0 = n mean(a_i b_i) and the conditional null
  variance V = n mean(v_i) (delta method, mobilink.nulls.analytic). Under an alternative with O/E = r, O is taken as
  approximately normal with mean r E0 and variance r V. The minimum detectable O/E at level alpha (two-sided) and the
  given power solves r E0 - z_power sqrt(r V) = E0 + z_alpha/2 sqrt(V) + 1/2 (enrichment), or the mirror image for
  depletion.
Net linked fraction at window w:
  An A gene of length g at a uniform position on a molecule of length L has a side that reaches w with probability
  (L - g - w)+ / (L - g), and the molecule holds lambda_A (L - g) such genes on average, so the expected number of A
  sides with room to reach w is 2 lambda_A sum_i (L_i - g - w)+ (averaged over the gene-length distribution). A side
  that meets a B gene before w ends there, so the number still at risk at w is that count times (1 - F_A(w)).
  Half-width (Greenwood-type, independent sides): 1.96 sqrt(F_A (1 - F_A) / n_w) / (1 - F_rnd), with n_w at risk at
  w, F_rnd = 1 - exp(-lambda_B w) and F_A = F_rnd + net (1 - F_rnd). It uses the smallest risk set, so it is mostly
  conservative: on the 8 G1b runs (bootstrap / formula half-width)^2 was 0.37-1.18 (median 0.7; step 1.14). deff
  scales the variance for dependence between sides (reads sharing a locus).
Depth enters through the number of molecules: n = bases / mean length.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import optimize, stats

from .nulls.analytic import conditional_variance, exposure, presence


def _sample_terms(lengths, len_a, len_b, lam_a, lam_b, min_frac):
    xa, xb = exposure(lengths, len_a, min_frac), exposure(lengths, len_b, min_frac)
    a, b = presence(xa, lam_a), presence(xb, lam_b)
    return float(np.mean(a * b)), conditional_variance(xa, xb, lam_a, lam_b) / len(lengths)


def detectable_oe(lengths, len_a, len_b, lam_a, lam_b, n_molecules, alpha=0.05, power=0.8, min_frac=0.8):
    """Minimum detectable O/E for enrichment and maximum detectable O/E for depletion (two-sided test)."""
    e1, v1 = _sample_terms(lengths, len_a, len_b, lam_a, lam_b, min_frac)
    e0, v = e1 * n_molecules, v1 * n_molecules
    za, zb = stats.norm.ppf(1 - alpha / 2), stats.norm.ppf(power)
    hi_crit = e0 + za * np.sqrt(v) + 0.5
    x = (zb * np.sqrt(v) + np.sqrt(zb**2 * v + 4 * e0 * hi_crit)) / (2 * e0)
    lo_crit = e0 - za * np.sqrt(v) - 0.5
    if lo_crit > 0:
        y = (-zb * np.sqrt(v) + np.sqrt(zb**2 * v + 4 * e0 * lo_crit)) / (2 * e0)
        r_lo = float(y**2)
    else:
        r_lo = np.nan  # depletion cannot be detected: even O = 0 is not significant
    return dict(expected=e0, var_null=v, oe_enrichment=float(x**2), oe_depletion=r_lo)


def power_at(lengths, len_a, len_b, lam_a, lam_b, n_molecules, oe, alpha=0.05, min_frac=0.8):
    """Approximate power of the two-sided test when the true O/E is `oe`."""
    e1, v1 = _sample_terms(lengths, len_a, len_b, lam_a, lam_b, min_frac)
    e0, v = e1 * n_molecules, v1 * n_molecules
    za = stats.norm.ppf(1 - alpha / 2)
    hi, lo = e0 + za * np.sqrt(v) + 0.5, e0 - za * np.sqrt(v) - 0.5
    sd = np.sqrt(oe * v)
    return float(stats.norm.sf((hi - oe * e0) / sd) + stats.norm.cdf((lo - oe * e0) / sd))


def sides_reaching(lengths, len_a, lam_a, w, n_molecules=None):
    """Expected number of A sides with room to reach w: 2 lambda_A sum_i (L_i - g - w)+, averaged over gene lengths g
    (whole genes). Scaled to n_molecules if given (lengths taken as a sample)."""
    L = np.asarray(lengths, float)
    s = 2 * lam_a * float(np.sum(exposure(np.clip(L - w, 0, None), len_a, 1.0)))
    return s if n_molecules is None else s * n_molecules / L.size


def chance_and_linked(net, lam_b, w):
    """(F_rnd, F_A) at w for a B call density lam_b (per bp) and a net linked fraction `net`."""
    f_r = 1 - np.exp(-lam_b * w)
    return float(f_r), float(f_r + net * (1 - f_r))


def at_risk(sides_with_room, net, lam_b, w):
    """Sides still at risk at w: those with room to reach w that have not met a B before w."""
    return sides_with_room * (1 - chance_and_linked(net, lam_b, w)[1])


def netlinked_halfwidth(n_sides_w, net, lam_b, w, deff=1.0):
    """Approximate 95% half-width of the net linked fraction at w from the number of A sides reaching w."""
    f_r, f_a = chance_and_linked(net, lam_b, w)
    return float(1.96 * np.sqrt(deff * f_a * (1 - f_a) / n_sides_w) / (1 - f_r))


def depth_table(
    lengths,
    len_a,
    len_b,
    lam_a,
    lam_b,
    lam_a_calls,
    lam_b_calls,
    bases,
    w=2000,
    net=0.1,
    deff=1.0,
    alpha=0.05,
    power=0.8,
    min_frac=0.8,
) -> pd.DataFrame:
    """Per depth (bases): molecules, expected co-occurrence, detectable O/E, A sides reaching w and the half-width of
    the net linked fraction.
    lam_a, lam_b: per-bp densities of the presence model (nulls.analytic.fit_density);
    lam_a_calls, lam_b_calls: calls per bp of exposure (calls_density), for the side count and the chance term."""
    L = np.asarray(lengths, float)
    rows = []
    for g in np.atleast_1d(bases):
        n = float(g) / L.mean()
        d = detectable_oe(L, len_a, len_b, lam_a, lam_b, n, alpha, power, min_frac)
        nw = at_risk(sides_reaching(L, len_a, lam_a_calls, w, n), net, lam_b_calls, w)
        rows.append(
            dict(
                bases=float(g),
                molecules=n,
                expected=d["expected"],
                oe_enrichment=d["oe_enrichment"],
                oe_depletion=d["oe_depletion"],
                sides_w=nw,
                halfwidth=netlinked_halfwidth(nw, net, lam_b_calls, w, deff) if nw > 0 else np.nan,
            )
        )
    return pd.DataFrame(rows)


def calls_density(n_calls: int, lengths, len_calls, min_frac=1.0) -> float:
    """Calls per bp of exposure: n_calls / sum_i x(L_i) (whole calls by default)."""
    return float(n_calls / np.sum(exposure(lengths, len_calls, min_frac)))


def required_bases(lengths, len_a, lam_a_calls, lam_b_calls, w, target_halfwidth, net=0.1, deff=1.0):
    """Bases needed for a net linked half-width at w of at most `target_halfwidth`."""
    L = np.asarray(lengths, float)

    def hw(g):
        n = g / L.mean()
        nw = at_risk(sides_reaching(L, len_a, lam_a_calls, w, n), net, lam_b_calls, w)
        return netlinked_halfwidth(nw, net, lam_b_calls, w, deff)

    return float(np.exp(optimize.brentq(lambda lg: hw(np.exp(lg)) - target_halfwidth, np.log(1e6), np.log(1e13))))


def depth_projection(test: pd.DataFrame, bases: float | None = None, targets=(0.05, 0.02, 0.01), alpha: float = 0.05,
                     power: float = 0.8) -> pd.DataFrame:
    """What a finished run can show for each pair, and what more sequence would give, from the test table of the run.
      halfwidth          half the width of the 95% interval of the test table, (net_w_hi - net_w_lo) / 2; when the
                         table has no interval, z(1 - alpha/2) * net_w_se;
      test_sd            |net_w / net_w_z|: the spread of the test statistic under chance, on the scale of the
                         estimate;
      detectable_net     (z(1 - alpha/2) + z(power)) * test_sd: the net linked fraction that the two-sided test at level
                         alpha finds with the given power at this depth (approximate; rough when the test is flagged
                         for few effective events);
      gbp_to_detect_<t>  bases * (detectable_net / t)^2 in Gbp, for each target t (only when `bases` is given): the
                         depth at which a net linked fraction of t becomes detectable, with the same read lengths and
                         community.
    The depth columns use the test, not the interval, on purpose. On eight public runs thinned to 1/2, 1/4 and 1/8 of
    their reads the variance of the test count fell in proportion to the depth (ratio 0.98 to 1.05, quartiles), and
    the projection from the thinned run to the full one was right within about 20% (quartiles 0.83 to 1.18). The
    delete-one-host interval did not follow: an 8-fold thinning widened it by a factor of 1.2 at the median, because
    it carries the differences between hosts, which more reads of the same sample do not remove. So more depth makes
    weaker linkage detectable but narrows the interval little. Where many reads already cover each locus they count
    once, and the depth columns are then the least that is needed."""
    za, zb = stats.norm.ppf(1 - alpha / 2), stats.norm.ppf(power)
    rows = []
    for _, r in test.iterrows():
        se, net, z = (float(r.get(k, np.nan)) for k in ("net_w_se", "net_w", "net_w_z"))
        sd = abs(net / z) if np.isfinite(net) and np.isfinite(z) and z != 0 and net != 0 else np.nan
        lo, hi = (float(r.get(k, np.nan)) for k in ("net_w_lo", "net_w_hi"))
        if not np.isfinite(net):
            hw = np.nan
        elif np.isfinite(lo) and np.isfinite(hi):  # the interval as reported: with few hosts it is wider than z * se
            hw = (hi - lo) / 2
        else:
            hw = za * se if np.isfinite(se) and se > 0 else np.nan
        row = dict(class_a=r["class_a"], class_b=r["class_b"], net_w=net, net_w_se=se, halfwidth=hw, test_sd=sd,
                   detectable_net=(za + zb) * sd)
        if bases is not None:
            row["gbp"] = bases / 1e9
            for t in targets:
                row[f"gbp_to_detect_{t:g}"] = bases / 1e9 * (row["detectable_net"] / t) ** 2
        rows.append(row)
    return pd.DataFrame(rows)
