"""Deterministic-equivalent (DE) formula sheet implementations for SCF Phase 1.

Conventions (see docs/model_card.md):
  - sigma_u^2 = 1 default; spike strengths l_j are eigenvalues of Lambda Lambda' / sigma_u^2.
  - Population eigenvalues of Sigma_X / sigma_u^2 are tau_j = 1 + l_j.
  - c = p/n.
  - gamma is given in factor coordinates (gamma_j multiplies population eigenvector u_j).

Every function here is pure numpy/scipy and carries its derivation in the
docstring. Each transcribed classical formula has a numerical self-check in
tests/test_identities.py before it is used anywhere (research plan WP 1.3).
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Classical random matrix facts (transcribed; self-checked in tests)
# ---------------------------------------------------------------------------


def mp_edges(c: float, sigma2: float = 1.0) -> tuple[float, float]:
    """Marchenko-Pastur bulk edges for sample covariance of white data.

    F6. For X (n x p) with iid columns of variance sigma2, the eigenvalues of
    X'X/n lie in [sigma2 (1-sqrt(c))^2, sigma2 (1+sqrt(c))^2] asymptotically,
    c = p/n. Source: Marchenko-Pastur (1967); Bai-Silverstein (2010) Ch. 3.
    """
    s = np.sqrt(sigma2)
    return (s * (1 - np.sqrt(c)) ** 2, s * (1 + np.sqrt(c)) ** 2)


def bbp_location(l: float, c: float, sigma2: float = 1.0) -> float:
    """BBP sample-eigenvalue location for a population spike of strength l.

    F4. Population eigenvalue tau = sigma2 (1 + l). For l > sqrt(c) the
    corresponding sample eigenvalue converges a.s. to

        mu(l) = tau (1 + c sigma2 / (tau - sigma2))
              = sigma2 (1 + l) (l + c) / l,

    and to the bulk edge (1+sqrt(c))^2 otherwise.

    CORRECTION vs research plan Section 2.3 item 4, which wrote
    mu(l) = (1+l)(1+cl)/l. The two agree only when l = 1 or c = 1
    (since 1 + cl = l + c iff (l-1)(c-1) = 0). The correct BBP form is
    (1+l)(1 + c/l); verified numerically in tests (test_bbp_location).
    Source: Baik-Ben Arous-Peche (2005), Thm 2.1 (real case); Johnstone (2001).
    """
    if l <= np.sqrt(c):
        return sigma2 * (1 + np.sqrt(c)) ** 2
    return sigma2 * (1 + l) * (l + c) / l


def bgn_overlap(l: float, c: float) -> float:
    """Squared eigenvector overlap |<u_j, v_hat_j>|^2 for a spike of strength l.

    F5. For l > sqrt(c): xi(l, c) = (1 - c/l^2) / (1 + c/l); else 0.
    Source: Benaych-Georges-Nadakuditi (2011), Adv. Math. 227(1), real case;
    also (2012) for the rectangular/singular-vector version.
    """
    if l <= np.sqrt(c):
        return 0.0
    return (1 - c / l**2) / (1 + c / l)


def tw_mu_sigma(n: int, p: int) -> tuple[float, float]:
    """Johnstone center/scale for the largest eigenvalue of X'X (unscaled).

    F7. For white X (n x p), entries N(0, sigma2), the largest eigenvalue of
    X'X (NOT divided by n) satisfies

        (lam_max - mu_np) / sigma_np  ->  TW1,

    mu_np   = (sqrt(n-1) + sqrt(p))^2,
    sigma_np = (sqrt(n-1) + sqrt(p)) * (1/sqrt(n-1) + 1/sqrt(p))^{1/3},

    multiplied by sigma2. Source: Johnstone (2001), Ann. Statist. 29(2).
    """
    sn, sp = np.sqrt(n - 1.0), np.sqrt(float(p))
    mu = (sn + sp) ** 2
    sigma = (sn + sp) * (1.0 / sn + 1.0 / sp) ** (1.0 / 3.0)
    return mu, sigma


# Tracy-Widom order-1 quantiles (standard published tables, e.g. Johnstone
# 2001 Table 1 / Edelman; used only as thresholds, calibrated by MC in tests).
TW1_Q95 = 0.9793
TW1_Q99 = 2.0234


def tw_threshold(n: int, p: int, sigma2: float = 1.0, q: float = TW1_Q99) -> float:
    """Upper-q threshold for lambda_max of X'X/n (white null), cov scale.

    F7. threshold = sigma2 * (mu_np + q * sigma_np) / n.
    """
    mu, sigma = tw_mu_sigma(n, p)
    return sigma2 * (mu + q * sigma) / n


# ---------------------------------------------------------------------------
# Bias functionals (M1). gamma given in factor coordinates; Q columns are the
# population eigenvectors u_j.
# ---------------------------------------------------------------------------


def sigma_x_eigs(l: np.ndarray, sigma2: float = 1.0) -> np.ndarray:
    """Eigenvalues of Sigma_X = sigma2 I + Lambda Lambda' (F1)."""
    l = np.asarray(l, float)
    return sigma2 * (1.0 + l)


def ols_bias_vector(l: np.ndarray, gamma: np.ndarray, sigma2: float = 1.0) -> np.ndarray:
    """Population OLS bias b_ols = Sigma_X^{-1} Lambda gamma (F1).

    Identity: Cov(X, Y) = Sigma_X beta + Lambda gamma with
    Sigma_X = sigma2 I + Lambda Lambda'. Since Lambda gamma =
    sum_j sigma_u sqrt(l_j) gamma_j u_j,

        plim beta_OLS - beta = sum_j [sqrt(l_j / sigma2) / (1 + l_j)] gamma_j u_j.

    EXACT finite-n fact (c <= 1): writing w = Lambda f + eps = X a + zeta with
    a = Sigma_X^{-1} Lambda gamma and zeta independent of X (joint Gaussianity,
    A2/A3), E[beta_OLS] - beta = a exactly, because (X'X)^{-1} X' X = I and
    E[(X'X)^{-1} X' zeta] = E[(X'X)^{-1} X'] E[zeta] = 0. So the simulated mean
    bias should equal this vector up to MC error at any n with p <= n.
    """
    l = np.asarray(l, float)
    gamma = np.asarray(gamma, float)
    return np.sqrt(l / sigma2) / (1.0 + l) * gamma


def ridge_bias_vector(
    l: np.ndarray, gamma: np.ndarray, lam: float, sigma2: float = 1.0
) -> np.ndarray:
    """Population ridge bias (F2): (Sigma_X + lam I)^{-1} Lambda gamma.

    plim beta_ridge - beta
        = sum_j [sigma_u sqrt(l_j) / (sigma_u^2 (1 + l_j) + lam)] gamma_j u_j,
    i.e., component form sqrt(l_j sigma2) / (sigma2 (1 + l_j) + lam) * gamma_j.
    Finite-n deviation order 1/n (shrinkage fluctuation of (Sigmahat+lam)^{-1}).
    """
    l = np.asarray(l, float)
    gamma = np.asarray(gamma, float)
    return np.sqrt(l * sigma2) / (sigma2 * (1.0 + l) + lam) * gamma


def pca_trim_bias_vector(
    l: np.ndarray, gamma: np.ndarray, k: int, sigma2: float = 1.0
) -> np.ndarray:
    """Population PCA-k trim bias (F3).

    Regressing Y on the top-k population PCs gives
        beta_trim = P_k beta + V_k (V_k' Sigma V_k)^{-1} V_k' Lambda gamma,
    so the bias keeps only the retained directions:
        bias = sum_{j <= k} [sqrt(l_j sigma2)/(sigma2(1+l_j))] gamma_j u_j.
    Trimming removes the bias of dropped components entirely and leaves the
    retained components untouched; the cost is signal loss
    1 - ||P_k beta||^2 / ||beta||^2 (large under dense beta, A4a).
    """
    l = np.asarray(l, float)
    gamma = np.asarray(gamma, float)
    out = np.zeros_like(gamma)
    out[:k] = np.sqrt(l[:k] * sigma2) / (sigma2 * (1.0 + l[:k])) * gamma[:k]
    return out


def mp_stieltjes_inv(lam: float, a: float) -> float:
    """E[1/(t + lam)] under the Marchenko-Pastur law with parameter a <= 1.

    F8 helper. For the nonzero spectrum of an n x n Gram matrix built from p
    columns (aspect a = n/p <= 1), the Stieltjes transform at z = -lam solves
    lam m^2 + (lam + a) m - 1 = 0, giving

        m_inv(lam) = (sqrt((lam+a)^2 + 4 lam) - (lam + a)) / (2 lam)  > 0,

    and E[1/(t+lam)] = m_inv. Source: standard MP Stieltjes equation
    (Bai-Silverstein 2010, Ch. 3); consistency checked numerically in tests.
    """
    disc = np.sqrt((lam + a) ** 2 + 4.0 * lam)
    return (disc - (lam + a)) / (2.0 * lam)


def ridge_capture(l: np.ndarray, lam: float, c: float, sigma2: float = 1.0) -> np.ndarray:
    """Ridge capture coefficients on spike directions, c > 1 (F8, PROVED).

    Shifted-resolvent assembly (docs/theory_T1_capture_law.md Section 6,
    closed 2026-08-25): with m_bar := lam-side companion Stieltjes value,

        m_bar(lam) = positive root of  lam*m^2 + (lam + c - 1)*m - 1 = 0
                   = [-(lam+c-1) + sqrt((lam+c-1)^2 + 4 lam)] / (2 lam),

        cap_j(lam) = (1 + l_j) m_bar / (1 + (1 + l_j) m_bar).

    Directional mean: E[<beta_hat_ridge - beta, q> | Q,beta] =
    -(cap_j(lam) - 1)<beta,q_j> + cap_j(lam) sqrt(l_j)/(1+l_j) gamma_j,
    i.e., the SAME structure as minnorm_capture with cap -> cap(lam).

    Collapse checks: lam->0 gives m_bar = 1/(c-1) and cap_j(0) =
    (1+l_j)/(c+l_j) (the proved min-norm capture law); lam->inf gives
    cap -> 0 (beta_hat -> 0). Reconciliation falsifier (2026-08-25):
    matches simulation to max |err| 0.003 over 9 (l,c,lam) cells while the
    superseded xi-split form errs up to 0.081 with lambda-dependent drift;
    see ridge_capture_superseded and docs/theory_T1_capture_law.md.
    """
    l = np.asarray(l, float)
    mb = mp_stieltjes_bulk_nside(float(lam), float(c))
    return (1.0 + l) * mb / (1.0 + (1.0 + l) * mb)


def mp_stieltjes_bulk_nside(lam: float, c: float) -> float:
    """lim tr(G + n*lam I)^{-1} for bulk Wishart K ~ W_n(I, df ~ c*n).

    Positive root of lam*m^2 + (lam + c - 1)*m - 1 = 0; equals
    int (x+lam)^{-1} d g(x) where g is the law of K/n eigenvalues
    (nonzero branch). Continuity anchor m_bar(0+) = 1/(c-1).
    Rationalized form 2/[(lam+c-1)+sqrt((lam+c-1)^2+4 lam)] avoids
    catastrophic cancellation as lam -> 0.
    """
    if lam <= 0:
        return 1.0 / (c - 1.0)
    s = np.sqrt((lam + c - 1.0) ** 2 + 4.0 * lam)
    return 2.0 / (lam + c - 1.0 + s)


def ridge_capture_superseded(
    l: np.ndarray, lam: float, c: float, sigma2: float = 1.0
) -> np.ndarray:
    """SUPERSEDED xi-split ridge capture (kept for audit trail).

    cap_j(lam) = xi_j nu_j/(nu_j+lam) + (1-xi_j)(1/c)(1-lam*m_inv(lam,1/c)).
    FALSIFIED 2026-08-25 by the decisive reconciliation check: its lam->0
    limit is the rejected superseded min-norm guess xi+(1-xi)/c, and
    simulation shows lambda-dependent drift up to 0.081 absolute (vs 0.003
    for the proved form) across (l,c,lam) cells at n=350, 200 reps.
    """
    l = np.asarray(l, float)
    xi = np.array([bgn_overlap(li, c) for li in l])
    nu = np.array([bbp_location(li, c, sigma2) for li in l])
    a = 1.0 / c
    bulk = (1.0 - lam * mp_stieltjes_inv(lam, a)) / c
    return xi * nu / (nu + lam) + (1.0 - xi) * bulk


def minnorm_capture(l: np.ndarray, c: float) -> np.ndarray:
    """Diagonal capture coefficients of E[P_row] on spike directions, c > 1 (F8).

    PILOT-VALIDATED CONJECTURE (Phase 1, WP 1.5): the capture law

        cap_j = (1 + l_j) / (c + l_j)

    matched simulation to ~0.5-1 percent at c = 5 across supercritical
    (l = 3 sqrt(c)) and subcritical (l = 0.5 sqrt(c), l ~ 0) components and
    both n = 400 and n = 2000; see docs/de_formula_sheet.md F8. Boundary
    anchors: l -> 0 gives the uniform rowspace fraction 1/c; l -> inf gives 1
    (consistent estimation of spike-aligned components). For c <= 1 OLS is
    full rank and capture is identically 1.

    The earlier xi-based formula xi(l,c) + (1 - xi(l,c))/c predicted 0.607 at
    (l, c) = (3 sqrt(5), 5) against a measured 0.657-0.659 and is recorded in
    bgn_capture_superseded for comparison only.
    """
    l = np.asarray(l, float)
    return (1.0 + l) / (c + l)


def bgn_capture_superseded(l: np.ndarray, c: float) -> np.ndarray:
    """Superseded xi-based capture guess (kept for audit trail)."""
    l = np.asarray(l, float)
    xi = np.array([bgn_overlap(li, c) for li in l])
    return xi + (1.0 - xi) / c


def minnorm_bias_vector(
    l: np.ndarray,
    gamma: np.ndarray,
    c: float,
    beta_spike: np.ndarray | None = None,
    sigma2: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Min-norm OLS bias decomposition for c > 1 (F8).

    Returns (bias_confounding, bias_total_extra_terms) where
        bias_confounding = sum_j capture_j * [sigma2 l_j/(sigma2(1+l_j))] gamma_j u_j,
        total bias adds the fit artifact
            sum_j (capture_j - 1) beta_j u_j + (1/c - 1) beta_perp
    with beta_spike the coordinates of beta on (u_j) (dense A4a beta:
    beta_j ~ N(0, 1/p), beta_perp norm^2 ~ 1 - r/p).

    Derivation: see minnorm_capture docstring; the zeta decomposition makes
    the confounding part exact at the DE level (only xi(l,c) itself is an
    asymptotic object). Validated against simulation in the pilot (c = 5
    overlay); flagged PROVISIONAL-NEW derivation in de_formula_sheet.md.
    """
    l = np.asarray(l, float)
    gamma = np.asarray(gamma, float)
    cap = minnorm_capture(l, c)
    bias_conf = cap * ols_bias_vector(l, gamma, sigma2)
    if beta_spike is None:
        return bias_conf, np.zeros_like(bias_conf)
    beta_spike = np.asarray(beta_spike, float)
    fit_spike = (cap - 1.0) * beta_spike
    return bias_conf, fit_spike


def minnorm_total_bias_norm(
    l: np.ndarray,
    gamma: np.ndarray,
    c: float,
    p: int,
    sigma2: float = 1.0,
) -> float:
    """RMS-over-beta prediction of ||E[beta_hat] - beta|| for c > 1 (F8).

    Under A4a, beta_j ~ N(0, 1/p)-scale coordinates and
    ||beta_perp||^2 ~ 1 - r/p, so

        E||bias||^2 = sum_j [(cap_j - 1)^2 / p + (cap_j * ols_j * gamma_j)^2]
                      + (1/c - 1)^2 (1 - r/p).
    """
    l = np.asarray(l, float)
    gamma = np.asarray(gamma, float)
    r = len(l)
    cap = minnorm_capture(l, c)
    ols = ols_bias_vector(l, gamma, sigma2)
    e2 = np.sum((cap - 1.0) ** 2 / p + (cap * ols * gamma) ** 2)
    e2 += (1.0 / c - 1.0) ** 2 * max(0.0, 1.0 - r / p)
    return float(np.sqrt(e2))


# ---------------------------------------------------------------------------
# Detection quantities (C2 scaffolding; frontier calibration is Phase 2)
# ---------------------------------------------------------------------------


def eff_detect_spike(
    l: np.ndarray, gamma: np.ndarray, c: float, sigma2: float = 1.0
) -> float:
    """Effective detection spike s_eff (F9, working definition).

    s_eff = || P_spike Sigma_X^{-1/2} Lambda gamma ||^2
          = sum_{j: l_j > sqrt(c)} [sigma2 l_j / (sigma2 (1 + l_j))] gamma_j^2.

    Only the supercritical-aligned part of b = Lambda gamma produces a coherent
    rank-r mean shift in the whitened cross-moment statistic S1; subcritical
    components hide in the bulk (s_eff = 0 at leading order). This is the
    formal content of "invisible" in the decoupling claim. The mapping from
    s_eff to power (the frontier curve) is calibrated in Phase 2 (OMH
    template); treated as heuristic in Phase 1.
    """
    l = np.asarray(l, float)
    gamma = np.asarray(gamma, float)
    sup = l > np.sqrt(c)
    return float(np.sum((sigma2 * l[sup] / (sigma2 * (1.0 + l[sup]))) * gamma[sup] ** 2))


# ---------------------------------------------------------------------------
# Phase 2 detection formulas: F12 (max-z law) and F13 (augmented BBP).
# Frozen spec: docs/detection_statistics.md; preregistration:
# docs/phase2_preregistration.md. NEW-DE items validated by WP 2.1 null cells
# before powering WP 2.3 gates.
# ---------------------------------------------------------------------------


def maxz_null_var(
    d: float,
    c: float,
    sigma_eps2: float = 1.0,
    sigma_y2: float | None = None,
) -> float:
    """Null variance of the standardized spike coordinate z_j (F12, ERRATUM 1).

    With Y standardized, b = X_c' Ytilde / n, w_j = v_j' b and
    z_j = sqrt(n) w_j / sqrt(d_j), under H0 (gamma = 0), A4a, A2/A3:

        Var(w_j) = d_j^2 / p + sigma_eps2 d_j / (n sigma_y2^2),
        =>  Var(z_j) = [ d_j / c + sigma_eps2 ] / sigma_y2^2,

    where sigma_y2 = beta' Sigma beta + sigma_eps2 (approximately
    tr(Sigma)/p + sigma_eps2 under A4a) enters because Y is standardized.
    ERRATUM: an earlier draft claimed Var(z_j) = 1 + c d_j; that inverted the
    rowspace factor (n d_j / p = d_j / c, not c d_j). Caught by
    tests/test_phase2.py::test_maxz_null_variance BEFORE any Phase 2 sweep
    data was generated (2026-08-24); see docs/detection_statistics.md
    Erratum 1. Practical reading: dense-beta leakage through spike directions
    GROWS like d_j/c, so calibrated spike-coordinate tests pay a large tax at
    small c; the empirical detection frontier must be read jointly with S1.
    """
    base = d / c + sigma_eps2
    return base / sigma_y2 ** 2 if sigma_y2 else base


def mp_quantile(q: float, c: float) -> float:
    """q-th quantile of the Marchenko-Pastur law (ratio c <= 1 branch).

    Density f(x) = sqrt((b - x)(x - a)) / (2 pi c x) on [a, b] with
    a = (1 - sqrt(c))^2, b = (1 + sqrt(c))^2. Bisection on the explicit CDF;
    used to convert a bulk-eigenvalue median into a sigma_u^2 estimate.
    """
    a, b = (1 - np.sqrt(c)) ** 2, (1 + np.sqrt(c)) ** 2

    def cdf(x):
        # numeric CDF via fine trapezoid on the explicit density
        xs = np.linspace(a, x, 4001)
        dens = np.sqrt(np.maximum((b - xs) * (xs - a), 0.0)) / (2 * np.pi * c * xs)
        return float(np.trapezoid(dens, xs))

    lo, hi = a, b
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if cdf(mid) < q:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def maxz_threshold(
    ktop: int, alpha: float = 0.05, two_sided: bool = True
) -> float:
    """Bonferroni Gaussian threshold for the max over ktop calibrated z's."""
    from scipy.stats import norm

    tail = alpha / (2 * ktop) if two_sided else alpha / ktop
    return float(norm.ppf(1.0 - tail))


def bbp_invert(nu: float, c: float) -> float:
    """Population spike l with sample location nu (F4 inverse).

    Solves (1 + l)(l + c)/l = nu for l > sqrt(c): l^2 + (1 + c - nu) l + c = 0,
    larger root. Returns np.nan if nu <= bulk edge (no supercritical spike can
    sit there). Used by the SEB tuner to map sample spikes back to l_j.
    """
    if nu <= (1 + np.sqrt(c)) ** 2:
        return np.nan
    b = 1.0 + c - nu
    disc = b * b - 4.0 * c
    if disc < 0:
        return np.nan
    return 0.5 * (-b + np.sqrt(disc))


def estimate_noise_scales(d_desc: np.ndarray, c: float):
    """Frozen observable-only estimates (se2_hat, sigma_y2_hat).

    se2_hat: sigma_u^2 proxy from the non-spike bulk median matched against
    the MP median (bulk median = sigma_u2 * MP_median(min(c,1))); capped at
    the observed bulk median. sigma_y2_hat = mean(d) + se2_hat
    (tr(Sigmahat)/p + sigma_eps2 under A4a). Used identically by the
    detection statistics and the SEB tuner so their calibrations agree.
    """
    d_desc = np.asarray(d_desc, float)
    edge = (1 + np.sqrt(c)) ** 2 * 1.05 if c <= 1 else np.inf
    tail = d_desc[d_desc < edge] if c <= 1 else d_desc
    med_bulk = float(np.median(tail)) if len(tail) else float(np.median(d_desc))
    se2_hat = med_bulk / mp_quantile(0.5, min(c, 1.0))
    se2_hat = float(min(se2_hat, med_bulk))
    sigma_y2 = float(np.mean(d_desc)) + se2_hat
    return se2_hat, sigma_y2


def aug_secular_term(lam: float, eigs_desc: np.ndarray) -> float:
    """sum_j tau_j^2 / (lam - tau_j) evaluated above the top eigenvalue.

    F13 helper: the resolvent trace entering the augmented secular equation.
    Caller guarantees lam > max(eigs); poles below lam are impossible then.
    """
    t = np.asarray(eigs_desc, float)
    return float(np.sum(t ** 2 / (lam - t)))


def aug_null_root(
    eigs_desc: np.ndarray, sigma_y2: float, p_eff: int | None = None
) -> float:
    """H0 location of lambda_max of the augmented moment matrix (F13).

    Under H0 the cross-moment vector b = Sigma beta / sigma_y has isotropic
    geometry (A4a), so E[b'(lam I - Sigma)^{-1} b]
      = sum_j tau_j^2 / (p sigma_y^2 (lam - tau_j)),
    and the out-of-bulk root solves

        lam = 1 + aug_secular_term(lam) / (p sigma_y2).

    Plug-in version uses estimated eigenvalues (descending, covariance scale)
    and sigma_y2_hat; see detection.py for the estimator. Bisection on
    [t_max + eps, upper]; the upper bracket comes from the conservative bound
    lam < 1 + S/(lam - t_max) solved as a quadratic.
    """
    t = np.sort(np.asarray(eigs_desc, float))[::-1]
    p = len(t)
    tmax = t[0]
    lo = tmax * (1.0 + 1e-9) + 1e-12

    def f(lam):
        return lam - 1.0 - aug_secular_term(lam, t) / (p * sigma_y2)

    # conservative upper bracket: replace every tau by tmax in denominators
    s_total = float(np.sum(t ** 2))
    # lam - 1 = S/(p sy2 (lam - tmax)) -> lam^2 - (1+tmax) lam + (tmax + S/(p sy2))
    disc = (1.0 + tmax) ** 2 - 4.0 * (tmax + s_total / (p * sigma_y2))
    hi = 0.5 * ((1.0 + tmax) + np.sqrt(max(disc, 1e-300))) * (1.0 + 1e-6) + 1e-6
    flo, fhi = f(lo), f(hi)
    if flo > 0:  # root coincides with the spike itself (degenerate); clamp
        return lo
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
        if hi - lo < 1e-10 * max(1.0, hi):
            break
    return 0.5 * (lo + hi)


def aug_h1_root(
    l: np.ndarray,
    gamma: np.ndarray,
    sigma_y2: float,
    sigma_u2: float = 1.0,
) -> float:
    """H1 location prediction of the augmented statistic (F13).

    Population-level secular equation with confounding:

        lam = 1 + [ sum_j tau_j^2 / (p sigma_y2 (lam - tau_j)) ]
                + [ sigma_u2 sum_j l_j gamma_j^2 / (sigma_y2 (lam - tau_j)) ],

    derived from E[b] = (Sigma beta + Lambda gamma)/sigma_y: the A4a beta part
    spreads isotropically while the gamma part contributes the coherent
    rank-r term. The O(1/sqrt(p)) cross term between beta and gamma is
    dropped (documented approximation, checked in WP 2.1 overlays).
    """
    l = np.asarray(l, float)
    gamma = np.asarray(gamma, float)
    tau = sigma_u2 * (1.0 + l)
    p = len(tau)

    def f(lam):
        beta_part = aug_secular_term(lam, tau) / (p * sigma_y2)
        conf_part = sigma_u2 * float(
            np.sum(l * gamma ** 2 / (lam - tau))
        ) / sigma_y2
        return lam - 1.0 - beta_part - conf_part

    tmax = float(tau.max())
    lo = tmax * (1.0 + 1e-9) + 1e-12

    def g(lam):
        return lam - 1.0 - aug_secular_term(lam, tau) / (p * sigma_y2)

    s_total = float(np.sum(tau ** 2)) + p * sigma_u2 * float(np.sum(l * gamma ** 2))
    disc = (1.0 + tmax) ** 2 - 4.0 * (tmax + s_total / (p * sigma_y2))
    hi = 0.5 * ((1.0 + tmax) + np.sqrt(max(disc, 1e-300))) * (1.0 + 1e-6) + 1e-6
    if f(lo) > 0 or g(hi) <= 0:
        return lo
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
        if hi - lo < 1e-10 * max(1.0, hi):
            break
    return 0.5 * (lo + hi)


def tw_width_cov_scale(n: int, p: int, sigma2: float = 1.0) -> float:
    """Johnstone TW width on the covariance scale for an (n, p) white design.

    Used as the v1 fluctuation scale of T_aug around its secular root with
    dims (n, p+1). Known to be only approximately right for spiked
    deformations; the MC-calibrated threshold variant is co-recorded per
    docs/detection_statistics.md.
    """
    _, sig = tw_mu_sigma(n, p)
    return sigma2 * sig / n


# ---------------------------------------------------------------------------
# Factor-number selection rules (baselines; constants flagged approximate)
# ---------------------------------------------------------------------------

# Onatski (2010) ED-ratio asymptotic critical values, k = 1..10 (his Table 1,
# standard normal quantile-based construction). APPROXIMATE transcription;
# used only for PCA-k baseline selection, refined in Phase 2.
ONATSKI_CRIT = np.array([2.19, 2.09, 2.04, 2.01, 1.99, 1.97, 1.96, 1.95, 1.94, 1.93])


def onatski_select(eigs_desc: np.ndarray, kmax: int = 10) -> int:
    """Onatski ratio rule: largest k with eigs[k-1]/eigs[k] > crit[k-1] (F10)."""
    eigs_desc = np.asarray(eigs_desc, float)
    kmax = min(kmax, len(eigs_desc) - 1, len(ONATSKI_CRIT))
    k = 0
    for j in range(kmax):
        if eigs_desc[j] <= 0 or eigs_desc[j + 1] <= 0:
            break
        if eigs_desc[j] / eigs_desc[j + 1] > ONATSKI_CRIT[j]:
            k = j + 1
        else:
            break
    return max(k, 0)


def bai_ng_select(eigs_desc: np.ndarray, n: int, p: int, kmax: int = 10) -> int:
    """Bai-Ng PC-style selector (F11). PROVISIONAL approximation.

    khat = argmin_k [ sum_{j>k} eigs_j + k * mean(eigs) * log(max(n,p)) * (n+p)/(np) ].
    The exact Bai-Ng (2002) penalty constants are not load-bearing for Phase 1
    (baseline selection only); flagged for exact transcription before Phase 2.
    """
    eigs_desc = np.asarray(eigs_desc, float)
    kmax = min(kmax, len(eigs_desc) - 1)
    scale = float(np.mean(eigs_desc)) * np.log(max(n, p)) * (n + p) / (n * p)
    vals = [np.sum(eigs_desc[k:]) + k * scale for k in range(0, kmax + 1)]
    return int(np.argmin(vals))


# ---------------------------------------------------------------------------
# Ledger hashing (mechanical verification hook, WP 1.1)
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Spectral weight families and SEB tuner targets (WP 2.2; frozen in
# docs/phase2_preregistration.md)
# ---------------------------------------------------------------------------


def trim_weights(d: np.ndarray, tau: float) -> np.ndarray:
    """Cevid et al. Trim transform on covariance eigenvalues: w = min(1, tau/d).

    Their eq. (3.3) caps SINGULAR values at tau; on the covariance scale
    d_cov = d_svd^2/n this is w_j = min(1, (tau_svd/svd_j)^2). We parameterize
    directly on the covariance scale with tau_cov; default tuning
    tau = median(d) corresponds exactly to their tau = median singular value.
    """
    d = np.asarray(d, float)
    return np.minimum(1.0, tau / np.maximum(d, 1e-300))


def lava_weights(d: np.ndarray, rho: float) -> np.ndarray:
    """LAVA / SDBoost spectral-loss weights w_j = 1/(1 + rho d_j).

    Nava et al. eq. (14): w_i = n lam2/(n lam2 + s_i^2) with s_i singular
    values; on the covariance scale s_i^2 = n d_i so w_i = 1/(1 + d_i/rho')
    with n lam2 = rho'. We absorb constants into rho >= 0.
    """
    d = np.asarray(d, float)
    return 1.0 / (1.0 + rho * d)


def sdboost_path_coefficients(
    z: np.ndarray, d: np.ndarray, w: np.ndarray, m: int, nu: float
) -> np.ndarray:
    """SDBoost linear-base-learner path alpha_j(m) = (z_j/d_j)(1-(1-nu w_j)^m).

    Direct consequence of their boosting recursion (Nava et al. Section 3):
    coordinate j of the coefficient vector in the right-singular basis moves
    toward its OLS value z_j/d_j at per-direction rate nu w_j. m -> inf
    recovers min-norm OLS; early stopping is what creates deconfounding.
    """
    d = np.maximum(np.asarray(d, float), 1e-300)
    return (np.asarray(z, float) / d) * (
        1.0 - (1.0 - nu * np.asarray(w, float)) ** int(m)
    )


def seb_predicted_mse(
    l_hat: np.ndarray,
    g2_hat: np.ndarray,
    d: np.ndarray,
    tau: float,
    c: float,
    n: int,
    sigma_eps2: float,
    sigma_u2: float = 1.0,
) -> float:
    """DE-predicted causal MSE of the soft-trim estimator at threshold tau.

    SEB tuner objective (preregistration WP 2.2): for weights
    w_j = min(1, tau/d_j),

        MSE_pred(tau) = sum_j [ w_j cap_j sqrt(l_j sigma2)/(sigma2(1+l_j)) ghat_j ]^2
                        + sum_j w_j^2 sigma_eps2 / (n d_j)
                        + sum_j (1 - w_j)^2 / p          (A4a signal loss)
                        + artifact floor (c > 1 rowspace proxy),

    where cap_j = minnorm_capture(l_j, c) at c > 1 and 1 otherwise, and
    l_hat/g2_hat are plug-in estimates from the spectrum and cross-moment
    coordinates. Pure function of observables when fed estimates; the ORACLE
    variant feeds the true (l, gamma^2) into the SAME objective (frozen
    definition of eb_oracle_tau, preregistration).
    """
    l = np.asarray(l_hat, float)
    g2 = np.asarray(g2_hat, float)
    d = np.asarray(d, float)
    w = trim_weights(d[: len(l)], tau)
    cap = minnorm_capture(l, c) if c > 1 else np.ones_like(l)
    coef = np.sqrt(l * sigma_u2) / (sigma_u2 * (1.0 + l))
    bias2 = float(np.sum((w * cap * coef * np.sqrt(g2)) ** 2))
    var_term = float(
        np.sum((w ** 2) * sigma_eps2 / (n * np.maximum(d[: len(l)], 1e-12)))
    )
    # A4a signal-attenuation cost: dense beta has ~1/p mass per direction,
    # systematic shrinkage (1 - w_j) contributes sum_j (1-w_j)^2 / p.
    atten = float(np.sum((1.0 - w) ** 2)) / max(c * n, 1)
    artifact = max(0.0, 1.0 / c - 1.0) ** 2 * float(np.mean(w)) if c > 1 else 0.0
    return bias2 + var_term + atten + artifact


def ucm_strength(eigs_desc: np.ndarray, n: int, p: int) -> float:
    """Rendsburg-et-al.-style confounding strength proxy (B2 baseline).

    Documented approximation of the UCM point estimate: the fraction of total
    variance carried by super-BBP-outlier directions relative to the white
    expectation. Uses only the spectrum; bootstrap thresholding is applied by
    the caller. NOT a faithful reimplementation of their full PE algorithm;
    flagged APPROXIMATE in docs/detection_statistics.md.
    """
    eigs_desc = np.asarray(eigs_desc, float)
    edge = (1.0 + np.sqrt(p / n)) ** 2
    excess = np.sum(np.maximum(eigs_desc - edge, 0.0))
    return float(excess / max(np.sum(eigs_desc), 1e-12))


def ledger_hash(doc_dir: str | Path) -> str:
    """sha256(model_card.md || assumption_ledger.md), first 12 hex chars."""
    d = Path(doc_dir)
    h = hashlib.sha256()
    for name in ("model_card.md", "assumption_ledger.md"):
        h.update((d / name).read_bytes())
    return h.hexdigest()[:12]
