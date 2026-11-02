"""Pure-python statistical engine for the scMR + CMap project.

No scipy dependency (scipy is not installed in this environment).
Implements: normal/t/chi2 distributions, IVW (FE + RE), MR-Egger with intercept,
weighted median, simple/weighted mode, Cochran Q, leave-one-out, F statistic,
Steiger directionality, coloc ABF (Giambartolomei 2014) and a simplified HEIDI.
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

# ----------------------------------------------------------------------------
# Special functions
# ----------------------------------------------------------------------------


def _betacf(a: float, b: float, x: float, max_iter: int = 300, eps: float = 3e-16) -> float:
    """Continued fraction for incomplete beta function (Numerical Recipes)."""
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def betainc(a: float, b: float, x: float) -> float:
    """Regularized incomplete beta I_x(a,b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    front = math.exp(math.log(x) * a + math.log(1.0 - x) * b - lbeta)
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def norm_cdf(z: float) -> float:
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def norm_sf(z: float) -> float:
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def norm_ppf(p: float) -> float:
    """Inverse normal CDF (Acklam's rational approximation)."""
    if p <= 0.0:
        return -math.inf
    if p >= 1.0:
        return math.inf
    a = [-3.969683028665376e01, 2.209460984245205e02, -2.759285104469687e02,
         1.383577518672690e02, -3.066479806614716e01, 2.506628277459239e00]
    b = [-5.447609879822406e01, 1.615858368580409e02, -1.556989798598866e02,
         6.680131188771972e01, -1.328068155288572e01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e00,
         -2.549732539343734e00, 4.374664141464968e00, 2.938163982698783e00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e00,
         3.754408661907416e00]
    plow, phigh = 0.02425, 1.0 - 0.02425
    if p < plow:
        q = math.sqrt(-2.0 * math.log(p))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    if p > phigh:
        q = math.sqrt(-2.0 * math.log(1.0 - p))
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
                ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    q = p - 0.5
    r = q * q
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / \
           (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)


def t_sf(t: float, df: float) -> float:
    """Upper tail of Student t."""
    if df <= 0:
        return float("nan")
    x = df / (df + t * t)
    p = 0.5 * betainc(df / 2.0, 0.5, x)
    return p if t > 0 else 1.0 - p


def t_two_sided(t: float, df: float) -> float:
    return min(1.0, 2.0 * t_sf(abs(t), df))


def chi2_sf(x: float, df: float) -> float:
    """Upper tail of chi-square via regularized upper incomplete gamma."""
    if x <= 0:
        return 1.0
    a, xx = df / 2.0, x / 2.0
    if xx < a + 1.0:
        # series for P(a,x)
        ap, s, delta = a, 1.0 / a, 1.0 / a
        for _ in range(1000):
            ap += 1.0
            delta *= xx / ap
            s += delta
            if abs(delta) < abs(s) * 1e-16:
                break
        p = s * math.exp(-xx + a * math.log(xx) - math.lgamma(a))
        return max(0.0, 1.0 - p)
    # continued fraction for Q(a,x)
    tiny = 1e-300
    b = xx + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b if abs(b) > tiny else 1.0 / tiny
    h = d
    for i in range(1, 1000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-16:
            break
    return min(1.0, math.exp(-xx + a * math.log(xx) - math.lgamma(a)) * h)


def fisher_combine(pvals: Sequence[float]) -> Tuple[float, float]:
    pv = [p for p in pvals if p is not None and 0 < p <= 1]
    if not pv:
        return float("nan"), float("nan")
    stat = -2.0 * sum(math.log(p) for p in pv)
    df = 2 * len(pv)
    return stat, chi2_sf(stat, df)


def bh_fdr(pvals: Sequence[float]) -> List[float]:
    """Benjamini-Hochberg adjusted p-values."""
    n = len(pvals)
    idx = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank in range(n, 0, -1):
        i = idx[rank - 1]
        val = min(prev, pvals[i] * n / rank)
        q[i] = val
        prev = val
    return q


# ----------------------------------------------------------------------------
# Linear algebra helpers (small, dependency free)
# ----------------------------------------------------------------------------


def wls_slope(x: Sequence[float], y: Sequence[float], w: Sequence[float]) -> Tuple[float, float, float]:
    """Weighted least squares y ~ x (no intercept). Returns (slope, se, intercept_free)."""
    sw = sum(w)
    swx2 = sum(wi * xi * xi for wi, xi in zip(w, x))
    swxy = sum(wi * xi * yi for wi, xi, yi in zip(w, x, y))
    if swx2 == 0:
        return float("nan"), float("nan"), float("nan")
    b = swxy / swx2
    se = math.sqrt(1.0 / swx2)
    return b, se, 0.0


def wls_line(x: Sequence[float], y: Sequence[float], w: Sequence[float]) -> Tuple[float, float, float, float]:
    """Weighted least squares with intercept. Returns (slope, slope_se, intercept, intercept_se)."""
    sw = sum(w)
    swx = sum(wi * xi for wi, xi in zip(w, x))
    swy = sum(wi * yi for wi, yi in zip(w, y))
    swxx = sum(wi * xi * xi for wi, xi in zip(w, x))
    swxy = sum(wi * xi * yi for wi, xi, yi in zip(w, x, y))
    det = sw * swxx - swx * swx
    if abs(det) < 1e-300:
        return float("nan"), float("nan"), float("nan"), float("nan")
    slope = (sw * swxy - swx * swy) / det
    intercept = (swxx * swy - swx * swxy) / det
    se_slope = math.sqrt(sw / det)
    se_int = math.sqrt(swxx / det)
    return slope, se_slope, intercept, se_int


# ----------------------------------------------------------------------------
# MR methods
# ----------------------------------------------------------------------------


def _clean(bx, by, sx, sy):
    out = []
    for a, b, c, d in zip(bx, by, sx, sy):
        if None in (a, b, c, d):
            continue
        if not all(math.isfinite(v) for v in (a, b, c, d)):
            continue
        if c <= 0 or d <= 0:
            continue
        out.append((a, b, c, d))
    return out


def ivw(bx, by, sx, sy, random_effects: bool = True) -> Dict[str, float]:
    d = _clean(bx, by, sx, sy)
    k = len(d)
    if k == 0:
        return {"nsnp": 0, "b": float("nan"), "se": float("nan"), "p": float("nan")}
    w = [1.0 / (sy[i] ** 2) for i in range(k)]
    b, se, _ = wls_slope([x[0] for x in d], [x[1] for x in d], w)
    if k > 1:
        q = sum(w[i] * (d[i][1] - b * d[i][0]) ** 2 for i in range(k))
        q_p = chi2_sf(q, k - 1)
        if random_effects and q > (k - 1) and q > 0:
            phi = q / (k - 1)
            se *= math.sqrt(phi)
    else:
        q, q_p = float("nan"), float("nan")
    z = b / se if se and math.isfinite(se) and se > 0 else float("nan")
    p = 2 * norm_sf(abs(z)) if math.isfinite(z) else float("nan")
    return {"nsnp": k, "b": b, "se": se, "p": p, "Q": q, "Q_pval": q_p, "z": z}


def mr_egger(bx, by, sx, sy) -> Dict[str, float]:
    d = _clean(bx, by, sx, sy)
    k = len(d)
    if k < 3:
        return {"nsnp": k, "b": float("nan"), "se": float("nan"), "p": float("nan")}
    bx_, by_, sx_, sy_ = [x[0] for x in d], [x[1] for x in d], [x[2] for x in d], [x[3] for x in d]
    # orient so that all bx > 0
    for i in range(k):
        if bx_[i] < 0:
            bx_[i] = -bx_[i]
            by_[i] = -by_[i]
    theta = [b / a for a, b in zip(bx_, by_)]
    theta.sort()
    med = theta[k // 2]
    sgn = 1.0 if med >= 0 else -1.0
    by2 = [y * sgn for y in by_]
    w = [1.0 / (sy_[i] ** 2 / bx_[i] ** 2) for i in range(k)]
    slope, se_s, intercept, se_i = wls_line(bx_, by2, w)
    df = k - 2
    t_int = intercept / se_i if se_i > 0 else float("nan")
    p_int = t_two_sided(t_int, df)
    # return slope with orientation restored
    b = slope * sgn
    se = se_s
    t_slope = slope / se_s if se_s > 0 else float("nan")
    p = t_two_sided(t_slope, df)
    return {"nsnp": k, "b": b, "se": se, "p": p,
            "egger_intercept": intercept * sgn, "egger_intercept_se": se_i,
            "egger_intercept_p": p_int}


def weighted_median(bx, by, sx, sy, n_boot: int = 1000, seed: int = 20261027) -> Dict[str, float]:
    import random
    d = _clean(bx, by, sx, sy)
    k = len(d)
    if k < 3:
        return {"nsnp": k, "b": float("nan"), "se": float("nan"), "p": float("nan")}
    bx_, by_, sx_, sy_ = [x[0] for x in d], [x[1] for x in d], [x[2] for x in d], [x[3] for x in d]
    for i in range(k):
        if bx_[i] < 0:
            bx_[i] = -bx_[i]
            by_[i] = -by_[i]
    ratio = [by_[i] / bx_[i] for i in range(k)]
    w = [bx_[i] ** 2 / sy_[i] ** 2 for i in range(k)]
    order = sorted(range(k), key=lambda i: ratio[i])
    ratio_s = [ratio[i] for i in order]
    w_s = [w[i] for i in order]
    tot = sum(w_s)
    cum, med = 0.0, ratio_s[0]
    for r, wi in zip(ratio_s, w_s):
        cum += wi
        if cum >= tot / 2.0:
            med = r
            break
    # bootstrap SE (parametric, first-order)
    rng = random.Random(seed)
    boots = []
    for _ in range(n_boot):
        br, bw = [], []
        for i in range(k):
            bxi = rng.gauss(bx_[i], sx_[i])
            byi = rng.gauss(by_[i], sy_[i])
            if bxi == 0:
                continue
            br.append(byi / bxi)
            bw.append(bxi ** 2 / sy_[i] ** 2)
        if not br:
            continue
        o = sorted(range(len(br)), key=lambda i: br[i])
        rs = [br[i] for i in o]
        ws = [bw[i] for i in o]
        t = sum(ws)
        c, m = 0.0, rs[0]
        for r, wi in zip(rs, ws):
            c += wi
            if c >= t / 2.0:
                m = r
                break
        boots.append(m)
    if len(boots) > 10:
        mu = sum(boots) / len(boots)
        se = math.sqrt(sum((x - mu) ** 2 for x in boots) / (len(boots) - 1))
    else:
        se = float("nan")
    z = med / se if se and se > 0 else float("nan")
    p = 2 * norm_sf(abs(z)) if math.isfinite(z) else float("nan")
    return {"nsnp": k, "b": med, "se": se, "p": p}


def _mode_estimator(bx, by, sy, phi: float = 1.0) -> float:
    """Hartwig 2017 weighted mode."""
    k = len(bx)
    ratio = [by[i] / bx[i] for i in range(k)]
    se_r = [abs(sy[i] / bx[i]) for i in range(k)]
    s = 1.0
    if k > 1:
        mu = sum(ratio) / k
        s = math.sqrt(sum((r - mu) ** 2 for r in ratio) / (k - 1))
        s = max(s, 1e-8)
    best_b, best_w = ratio[0], -1.0
    for j in range(k):
        w = 0.0
        for i in range(k):
            sd = se_r[i] * phi
            if sd <= 0:
                continue
            u = (ratio[i] - ratio[j]) / sd
            w += norm_cdf(1.0 / s * u) - norm_cdf(-1.0 / s * u) if s > 0 else 0.0
        if w > best_w:
            best_w, best_b = w, ratio[j]
    return best_b


def simple_mode(bx, by, sy) -> Dict[str, float]:
    d = _clean(bx, by, [1.0] * len(bx), sy)
    k = len(d)
    if k < 3:
        return {"nsnp": k, "b": float("nan"), "se": float("nan"), "p": float("nan")}
    bx_, by_, sy_ = [x[0] for x in d], [x[1] for x in d], [x[3] for x in d]
    b = _mode_estimator(bx_, by_, sy_, 1.0)
    # bootstrap
    import random
    rng = random.Random(11)
    boots = []
    for _ in range(500):
        idx = [rng.randrange(k) for _ in range(k)]
        boots.append(_mode_estimator([bx_[i] for i in idx], [by_[i] for i in idx], [sy_[i] for i in idx], 1.0))
    se = math.sqrt(sum((x - b) ** 2 for x in boots) / max(1, len(boots) - 1)) if len(boots) > 5 else float("nan")
    z = b / se if se and se > 0 else float("nan")
    p = 2 * norm_sf(abs(z)) if math.isfinite(z) else float("nan")
    return {"nsnp": k, "b": b, "se": se, "p": p}


def weighted_mode(bx, by, sx, sy) -> Dict[str, float]:
    d = _clean(bx, by, sx, sy)
    k = len(d)
    if k < 3:
        return {"nsnp": k, "b": float("nan"), "se": float("nan"), "p": float("nan")}
    bx_, by_, sy_ = [x[0] for x in d], [x[1] for x in d], [x[3] for x in d]
    ivw_res = ivw(bx_, by_, [x[2] for x in d], sy_, random_effects=True)
    phi = 1.0
    if ivw_res.get("Q") and k > 1 and ivw_res["Q"] > 0:
        phi = max(1.0, ivw_res["Q"] / (k - 1))
    b = _mode_estimator(bx_, by_, sy_, phi)
    import random
    rng = random.Random(13)
    boots = []
    for _ in range(500):
        idx = [rng.randrange(k) for _ in range(k)]
        boots.append(_mode_estimator([bx_[i] for i in idx], [by_[i] for i in idx], [sy_[i] for i in idx], phi))
    se = math.sqrt(sum((x - b) ** 2 for x in boots) / max(1, len(boots) - 1)) if len(boots) > 5 else float("nan")
    z = b / se if se and se > 0 else float("nan")
    p = 2 * norm_sf(abs(z)) if math.isfinite(z) else float("nan")
    return {"nsnp": k, "b": b, "se": se, "p": p}


def leave_one_out(bx, by, sx, sy) -> List[Dict[str, float]]:
    d = _clean(bx, by, sx, sy)
    k = len(d)
    out = []
    if k < 3:
        return out
    for i in range(k):
        sub = [d[j] for j in range(k) if j != i]
        r = ivw([x[0] for x in sub], [x[1] for x in sub], [x[2] for x in sub], [x[3] for x in sub])
        r["dropped"] = i
        out.append(r)
    return out


def f_statistic(bx, sx, n: int) -> float:
    if not bx:
        return float("nan")
    k = len(bx)
    if k == 0:
        return float("nan")
    r2 = sum((b * b) / (b * b + n * s * s) for b, s in zip(bx, sx)) / k
    if r2 <= 0 or r2 >= 1:
        return float("nan")
    return (n - k - 1) / k * r2 / (1 - r2)


def steiger(b_exposure: float, b_outcome: float, n_exposure: int, n_outcome: int,
            ncase: Optional[int] = None, ncontrol: Optional[int] = None) -> Dict[str, float]:
    """Approximate Steiger directionality test (r2 comparison)."""
    def r2_of(b, n, binary=False, ncase_=None, ncontrol_=None):
        if binary and ncase_ and ncontrol_:
            # convert log-OR to correlation via liability-scale approximation
            return b * b * 0.25 / (b * b * 0.25 + math.pi ** 2 / 3)
        denom = b * b + n
        return b * b / denom if denom else float("nan")
    r2e = r2_of(b_exposure, n_exposure)
    r2o = r2_of(b_outcome, n_outcome, binary=ncase is not None, ncase_=ncase, ncontrol_=ncontrol)
    return {"r2_exposure": r2e, "r2_outcome": r2o, "correct_direction": bool(r2e > r2o)}


# ----------------------------------------------------------------------------
# Colocalization (coloc.abf, Giambartolomei et al. 2014)
# ----------------------------------------------------------------------------

_DEFAULT_PRIORS = (1e-4, 1e-4, 1e-5)


def _log_abf(z: float, se: float, sd_prior: float = 0.15) -> float:
    """Wakefield log approximate Bayes factor for a single variant.

    lbf = 0.5 * ( log(1 - r) + r * z^2 ),  r = W / (W + V),
    with W = sd_prior^2 the prior variance of the effect and V = se^2 the
    squared standard error of the observed effect. Computed in log space so
    that strong signals (z ~ 15) do not overflow.
    """
    if z is None or se is None or not math.isfinite(z) or not math.isfinite(se):
        return 0.0
    if se <= 0:
        return 0.0
    w = sd_prior ** 2
    v = se ** 2
    r = w / (w + v)
    if r <= 0 or r >= 1:
        return 0.0
    return 0.5 * (math.log(1.0 - r) + r * z * z)


def _logsumexp(vals: Sequence[float]) -> float:
    vals = [v for v in vals if v is not None and math.isfinite(v)]
    if not vals:
        return float("-inf")
    m = max(vals)
    return m + math.log(sum(math.exp(v - m) for v in vals))


def coloc_abf(beta1, se1, beta2, se2, n1: Optional[int] = None, n2: Optional[int] = None,
              p12: float = 1e-5, sd_prior: float = 0.15) -> Dict[str, float]:
    """coloc.abf on two sets of summary statistics aligned on the same variants.

    Implements Giambartolomei et al. (2014) with the Wakefield approximate
    Bayes factor, evaluated in log space:

        l1 = sum_i log ABF_1(i),  l2 = sum_i log ABF_2(i),
        l12 = sum_i log ABF_1(i) + log ABF_2(i)

        H0 = 1
        H1 = p1 * exp(l1)
        H2 = p2 * exp(l2)
        H3 = p1 * p2 * exp(l1 + l2)
        H4 = p12 * exp(l12)

    Note the H3 term uses the product of the two marginal sums, not the
    per-variant product; only H4 sums per-variant joint ABFs.
    """
    d = []
    for b1, s1, b2, s2 in zip(beta1, se1, beta2, se2):
        if None in (b1, s1, b2, s2):
            continue
        if not all(math.isfinite(v) for v in (b1, s1, b2, s2)):
            continue
        if s1 <= 0 or s2 <= 0:
            continue
        d.append((b1, s1, b2, s2))
    n = len(d)
    if n < 2:
        return {"nsnp": n, "PP.H0": float("nan"), "PP.H1": float("nan"), "PP.H2": float("nan"),
                "PP.H3": float("nan"), "PP.H4": float("nan")}
    lb1 = [_log_abf(b1 / s1, s1, sd_prior) for b1, s1, _, _ in d]
    lb2 = [_log_abf(b2 / s2, s2, sd_prior) for _, _, b2, s2 in d]
    p1, p2, p12v = _DEFAULT_PRIORS
    l1 = _logsumexp(lb1)
    l2 = _logsumexp(lb2)
    l12 = _logsumexp([a + b for a, b in zip(lb1, lb2)])
    lh = [0.0, math.log(p1) + l1, math.log(p2) + l2,
          math.log(p1) + math.log(p2) + l1 + l2, math.log(p12v) + l12]
    norm = _logsumexp(lh)
    post = [math.exp(v - norm) for v in lh]
    return {"nsnp": n, "PP.H0": post[0], "PP.H1": post[1], "PP.H2": post[2],
            "PP.H3": post[3], "PP.H4": post[4]}


def heidi(b_eqtl: Sequence[float], b_gwas: Sequence[float], se_gwas: Sequence[float],
          top_index: int = 0, ld_r2: Optional[Sequence[float]] = None) -> Dict[str, float]:
    """Simplified HEIDI heterogeneity test across cis variants.

    Tests whether the eQTL-GWAS effect ratio is constant across SNPs
    (pleiotropy / two-distinct-causal-variants signal). Uses only SNPs in
    moderate LD with the top signal when ld_r2 is supplied.
    """
    pairs = []
    for i, (be, bg, sg) in enumerate(zip(b_eqtl, b_gwas, se_gwas)):
        if be is None or bg is None or sg is None:
            continue
        if not all(math.isfinite(v) for v in (be, bg, sg)):
            continue
        if abs(be) < 1e-8 or sg <= 0:
            continue
        if ld_r2 is not None and i < len(ld_r2) and ld_r2[i] is not None and ld_r2[i] < 0.05:
            continue
        pairs.append((bg / be, sg / abs(be)))
    if len(pairs) < 3:
        return {"nsnp": len(pairs), "Q": float("nan"), "heidi_p": float("nan")}
    d = [p[0] for p in pairs]
    s = [p[1] for p in pairs]
    w = [1.0 / (x ** 2) for x in s]
    sw = sum(w)
    bbar = sum(wi * di for wi, di in zip(w, d)) / sw
    q = sum(wi * (di - bbar) ** 2 for wi, di in zip(w, d))
    df = len(d) - 1
    return {"nsnp": len(d), "Q": q, "heidi_p": chi2_sf(q, df), "b_smr": bbar}


def smr_wald(b_eqtl: float, se_eqtl: float, b_gwas: float, se_gwas: float) -> Dict[str, float]:
    if not all(math.isfinite(v) for v in (b_eqtl, se_eqtl, b_gwas, se_gwas)):
        return {"b_smr": float("nan"), "se_smr": float("nan"), "p_smr": float("nan")}
    if abs(b_eqtl) < 1e-12:
        return {"b_smr": float("nan"), "se_smr": float("nan"), "p_smr": float("nan")}
    b = b_gwas / b_eqtl
    se = se_gwas / abs(b_eqtl)
    z = b / se if se > 0 else float("nan")
    p = 2 * norm_sf(abs(z)) if math.isfinite(z) else float("nan")
    return {"b_smr": b, "se_smr": se, "p_smr": p}


# ----------------------------------------------------------------------------
# CMap scoring helpers
# ----------------------------------------------------------------------------


def connectivity_score(up_overlap: Sequence[str], dn_overlap: Sequence[str],
                       n_up_query: int, n_dn_query: int,
                       n_up_ref: int = 150, n_dn_ref: int = 150) -> float:
    """Weighted connectivity score in the spirit of the CMap tau metric."""
    def ks(ov, nq, nr):
        if nq == 0 or nr == 0:
            return 0.0
        obs = len(ov)
        exp = nq * nr / 20000.0
        return (obs - exp) / math.sqrt(max(exp, 1e-9))
    a = ks(up_overlap, n_up_query, n_up_ref)
    b = ks(dn_overlap, n_dn_query, n_dn_ref)
    return a - b


def tau_score(scores: Sequence[float]) -> List[float]:
    """Normalize raw connectivity scores to the [-100, 100] tau scale."""
    out = []
    for s in scores:
        if s >= 0:
            mx = max(scores) if max(scores) > 0 else 1.0
            out.append(100.0 * s / mx)
        else:
            mn = min(scores) if min(scores) < 0 else -1.0
            out.append(-100.0 * s / mn)
    return out
