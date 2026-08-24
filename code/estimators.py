"""SCF Phase 2 estimator suite (WP 2.2). Frozen roster:
docs/phase2_preregistration.md. All estimators consume the shared spectrum
(descending eigenpairs of X_c'X_c/n) so per-rep cost is one eigendecomposition
plus cheap per-estimator algebra.

Every estimator returns a coefficient vector in R^p. Tuning budgets are fixed
module constants (logged into parquet rows by the runner) per the fair
comparison protocol (research plan Section 8.3).

Sources pinned for baselines:
- Cevid, Buhlmann, Meinshausen (AoS 2020), arXiv:1811.05352: Trim transform
  d~_i = min(d_i^svd, median(d^svd)), F = U diag(d~/d) U', regression on
  transformed data (their eqs. 3.1-3.3).
- Nava, Buhlmann, Sigrist (2026), arXiv:2607.09371: LAVA-type spectral loss
  w_i = sigma_e2/(sigma_r2 s_i^2 + sigma_e2); variance components by marginal
  likelihood ell(theta) (their Section 4.2.1); linear-base-learner boosting
  path alpha_j(m) = (z_j/d_j)(1-(1-nu w_j)^m) (their Section 3 recursion);
  stopping by BLUP-corrected K-fold CV (their Section 4.2.2).
"""
from __future__ import annotations

import math

import numpy as np

from de_formulas import (
    bbp_invert,
    lava_weights,
    seb_predicted_mse,
    sdboost_path_coefficients,
    trim_weights,
)

RIDGE_LAM_GRID = (0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0)
SDBOOST_NU = 0.1
SDBOOST_M_GRID = tuple(int(m) for m in np.unique(np.round(np.logspace(0, 3.2, 14))))
EB_TAU_FRAC_GRID = tuple(
    float(f) for f in np.concatenate([np.linspace(0.05, 0.95, 19), [1.5, 2.0, 4.0]])
)
N_FOLDS = 5


# ---------------------------------------------------------------------------
# shared pieces
# ---------------------------------------------------------------------------


def center_columns(X: np.ndarray) -> np.ndarray:
    return X - X.mean(axis=0, keepdims=True)


def standardize_response(Y: np.ndarray) -> np.ndarray:
    y = Y - Y.mean()
    sd = float(np.std(Y)) or 1.0
    return y / sd


def spectral_fit(d: np.ndarray, V: np.ndarray, Xty_over_n: np.ndarray, w: np.ndarray):
    """beta = V diag(w_j / d_j) V' (X'Y/n); w_j = weight on covariance scale."""
    coef = V.T @ Xty_over_n
    return V @ ((w * coef) / d)


def kfold_indices(n: int, k: int, rng: np.random.Generator):
    idx = rng.permutation(n)
    return np.array_split(idx, k)


# ---------------------------------------------------------------------------
# roster members; signature (Xc, Ys, eig, rng) -> beta_hat
# Xc column-centered design, Ys standardized response, eig = (d desc, V)
# ---------------------------------------------------------------------------


def est_ols(Xc, Ys, eig, rng=None):
    d, V = eig
    return fit_min_norm(Xc, Ys, eig)


def fit_min_norm(Xc, Ys, eig):
    d, V = eig
    return spectral_fit(d, V, Xc.T @ Ys / len(Ys), np.ones_like(d))


def loo_scores_spectral(
    U: np.ndarray, Ys: np.ndarray, d: np.ndarray, weight_fun, grid
):
    """Exact leave-one-out CV scores for any diagonal spectral fit.

    Any member of the family beta = V diag(w_j/d_j) U'Y induces the fitted
    operator Yhat = U diag(w_j) U' Y, whose hat-matrix diagonal is
    h_ii(w) = sum_j w_j u_ij^2. The exact LOO identity
    r_i^(-i) = (y_i - yhat_i) / (1 - h_ii(w)) then gives
    LOO-MSE(w) = mean(r_i^(-i)^2) at O(n k) per grid point, replacing
    fold-wise eigendecompositions (documented equivalence; budgets logged).
    Ridge corresponds to w_j = d_j/(d_j + lam) on the covariance scale.
    """
    n = len(Ys)
    z = U.T @ Ys
    U2 = U ** 2
    scores = []
    for param in grid:
        w = weight_fun(param)
        yhat = U @ (w * z)
        h = U2 @ w
        denom = np.maximum(1.0 - h, 1e-12)
        r_loo = (Ys - yhat) / denom
        scores.append(float(np.mean(r_loo ** 2)))
    return scores


def left_factors(Xc: np.ndarray, eig):
    """Left singular vectors U (n x k) from the shared spectrum."""
    d, V = eig
    return Xc @ V / np.sqrt(len(Xc) * d)[None, :]


def _cv_ridge(Xc, Ys, eig, rng, lam_grid=RIDGE_LAM_GRID):
    """Ridge with lambda by EXACT LOO-CV over the frozen grid (see
    loo_scores_spectral). Returns (beta_hat, info)."""
    n = len(Ys)
    d, V = eig
    U = left_factors(Xc, eig)

    def wf(lam):
        return d / (d + lam)

    scores = loo_scores_spectral(U, Ys, d, wf, lam_grid)
    lam = float(lam_grid[int(np.argmin(scores))])
    return _ridge_on_spectrum(d, V, Xc, Ys, lam), {"lam": lam, "grid": len(lam_grid)}


def _ridge_on_spectrum(d, V, Xc, Ys, lam):
    n = len(Ys)
    rhs = V.T @ (Xc.T @ Ys / n)
    return V @ (rhs / (d + lam))


def est_ridge_cv(Xc, Ys, eig, rng):
    return _cv_ridge(Xc, Ys, eig, rng)[0]


def est_pca_k(Xc, Ys, eig, rng=None, k: int = 0):
    d, V = eig
    w = np.zeros_like(d)
    w[:k] = 1.0
    return spectral_fit(d, V, Xc.T @ Ys / len(Ys), w)


def est_cevid_default(Xc, Ys, eig, rng=None):
    """Trim transform with tau = median singular value + OLS on transformed data.

    On the covariance scale their singular-value cap at median(d_svd) is
    w_j = min(1, median(d_svd)^2 / (n d_cov_j)); the global factor 1/n is
    absorbed by the fit invariance, so w_j = min(1, median_svd^2/(n d_j)).
    """
    n = len(Ys)
    d, V = eig
    med_svd2 = n * float(np.median(d))
    w = trim_weights(d * n, med_svd2)
    return spectral_fit(d, V, Xc.T @ Ys / n, w)


def est_lava_transform_ols(Xc, Ys, eig, rng=None, rho: float | None = None):
    """LAVA transform at rho = n*lam2 chosen by their lambda2 rule proxy.

    Their suggested rule sets n lam2 near bulk scale (sqrt(n)+sqrt(p))^2;
    we use exactly that as the default tuning (documented transcription of
    their Section 3 discussion via Nava et al. eq. 12 context).
    """
    n, p = Xc.shape
    d, V = eig
    if rho is None:
        rho = float((math.sqrt(n) + math.sqrt(p)) ** 2)
    w = lava_weights(d, rho)
    return spectral_fit(d, V, Xc.T @ Ys / n, w)


def sdboost_marginal_ll(theta_log: np.ndarray, s: np.ndarray, z: np.ndarray,
                        perp2: float, n: int) -> float:
    sr2, se2 = np.exp(theta_log)
    var = se2 + sr2 * s
    ll = -0.5 * np.sum(np.log(var) + z ** 2 / var)
    ll += -0.5 * (n - len(s)) * math.log(se2) - perp2 / (2 * se2)
    return ll


def fit_sdboost_linear_eb(Xc, Ys, eig, rng, nu: float = SDBOOST_NU):
    """SDBoost linear special case: EB(LAVA weights) + BLUP-corrected CV stop.

    Faithful composition of Nava et al.: (1) variance components by marginal
    likelihood on the XX' spectrum (their 4.2.1, single-shot version: the
    alternating theta/f updates are initialized at the ridge-small fit);
    (2) boosting path coefficients (their Section 3 recursion, exact);
    (3) stopping time by BLUP-corrected K-fold CV (their 4.2.2 formula).
    Returns (beta_hat, info dict).
    """
    from scipy.optimize import minimize

    if rng is None:
        rng = np.random.default_rng(0)
    n, p = Xc.shape
    d, V = eig
    U = Xc @ V / np.sqrt(n * d)[None, :]
    z = U.T @ Ys  # coordinates along left singular vectors (n-side)
    # variance components: initialize around ridge-small residual split
    lam_small = 0.01 * float(np.median(d))
    resid = Ys - _ridge_on_spectrum(d, V, Xc, Ys, lam_small) @ Xc.T
    perp2 = float(resid @ resid)
    s_full = n * d
    # coarse 2-D grid over log variance components, then local refine
    grid_sr = np.linspace(-16.0, 6.0, 23)
    grid_se = np.linspace(-16.0, 6.0, 23)
    ll_vals = np.array([
        [-sdboost_marginal_ll(np.array([a, b]), s_full, z, perp2, n)
         for b in grid_se] for a in grid_sr
    ])
    i0 = np.unravel_index(np.argmin(ll_vals), ll_vals.shape)
    x0 = np.array([grid_sr[i0[0]], grid_se[i0[1]]])
    from scipy.optimize import minimize

    res = minimize(
        lambda t: -sdboost_marginal_ll(t, s_full, z, perp2, n),
        x0, method="L-BFGS-B", bounds=[(-30.0, 10.0), (-30.0, 10.0)],
    )
    sr2, se2 = np.exp(res.x)

    def path_beta(m, d_svd, w, z_tr):
        alpha = sdboost_path_coefficients(z_tr, d_svd, w, m, nu)
        return V @ alpha

    # BLUP-corrected CV over the iteration grid. Frozen approximation: fold
    # random-effect operators use the FULL-data left singular subspace U
    # (rows masked to the training fold) instead of refolding spectra; the
    # kinship operator is dominated by the top-left-singular subspace, and
    # this removes all fold-wise eigendecompositions (documented in the
    # preregistration deviation log if it materially changes results).
    n_folds = N_FOLDS
    folds = kfold_indices(n, n_folds, rng)
    # cache fold-invariant quantities: cross Grams and training coordinates
    fold_cache = []
    for f in folds:
        mask = np.ones(n, bool)
        mask[f] = False
        Xtr = Xc[mask]
        n_tr = int(mask.sum())
        fold_cache.append({
            "mask": mask, "Xtr": Xtr, "n_tr": n_tr,
            "z_tr": U[mask].T @ Ys[mask],
            "cross": Xc[f] @ Xtr.T,
        })
    cv_scores = []
    for m in SDBOOST_M_GRID:
        tot = 0.0
        for fc in fold_cache:
            mask, Xtr, n_tr = fc["mask"], fc["Xtr"], fc["n_tr"]
            Ytr = Ys[mask]
            w_tr = 1.0 / (1.0 + (sr2 / se2) * (n_tr * d))
            bm = path_beta(m, np.sqrt(n_tr * d), w_tr, fc["z_tr"])
            rtr = Ytr - Xtr @ bm
            # Woodbury on Sig_tr = se2 I + sr2 * Utr S_tr Utr', shared-U approx
            s_tr = n_tr * d
            coef_u = U[mask].T @ rtr
            shrink = (sr2 * s_tr) / (se2 + sr2 * s_tr)
            sol = (rtr - U[mask] @ (shrink * coef_u)) / se2
            blup = sr2 * (fc["cross"] @ sol)
            pred = Xc[~mask] @ bm + blup
            tot += float(np.mean((Ys[~mask] - pred) ** 2))
        cv_scores.append(tot / n_folds)
    m_star = SDBOOST_M_GRID[int(np.argmin(cv_scores))]
    w_full = 1.0 / (1.0 + (sr2 / se2) * (n * d))
    beta = path_beta(m_star, np.sqrt(n * d), w_full, z)
    info = {"sr2": sr2, "se2": se2, "m": m_star, "grid": len(SDBOOST_M_GRID)}
    return beta, info


def spectrum_of(Xc: np.ndarray):
    """Descending eigenpairs of X_c'X_c/n reusing simulator.spectrum."""
    from simulator import spectrum

    return spectrum(Xc)


# ---------------------------------------------------------------------------
# OUR estimator: SEB-tuned soft trim
# ---------------------------------------------------------------------------


def estimate_spike_profile(d: np.ndarray, c: float, kmax: int = 10):
    """Map sample eigenvalues to population spike estimates l_hat_j.

    Walk down the top kmax eigenvalues; while bbp_invert yields l > sqrt(c),
    record it; stop at the first non-outlier. Subcritical directions get
    l_hat = small floor (1e-3) so capture ~ 1/c-ish behavior is preserved.
    """
    l_hats = []
    for dj in d[:kmax]:
        lj = bbp_invert(dj, c)
        if np.isnan(lj):
            break
        l_hats.append(lj)
    return np.array(l_hats if l_hats else [1e-3])


def estimate_gamma2(
    Xc: np.ndarray, Ys: np.ndarray, eig, l_hat: np.ndarray, c: float
) -> np.ndarray:
    """Cross-moment mixture estimate of gamma_j^2 (SEB plug-in).

    z_j = sqrt(n) v_j'b / sqrt(d_j) has null variance
    var_cal = (se2 + d_j/c)/sigma_y2 (F12 Erratum 1) and H1 mean
    mu_j = sqrt(n l_j) g dir_j / (sqrt(1+l_j) sigma_y). Positive-part moment
    inversion:

        ghat_j^2 = max(0, z_j^2 - var_cal) (1 + l_hat_j) sigma_y^2 / (n l_hat_j),

    conservative since dir^2 <= 1 and chi-square fluctuations inflate it.
    Frozen as the SEB estimator; scales from estimate_noise_scales so
    detection and tuning agree.
    """
    n = len(Ys)
    d, V = eig
    from de_formulas import estimate_noise_scales

    se2, sigma_y2 = estimate_noise_scales(d, c)
    sd_y = float(np.std(Ys)) or 1.0
    b = Xc.T @ Ys / (n * sd_y)
    out = []
    for j, lj in enumerate(l_hat):
        wj = float(V[:, j] @ b)
        zj2 = n * wj ** 2 / d[j]
        var_cal = (se2 + d[j] / c) / sigma_y2 ** 2
        excess = max(0.0, zj2 - var_cal)
        denom = n * max(lj, 1e-6)
        out.append(excess * (1.0 + lj) * sigma_y2 / denom)
    return np.array(out)


def _noise_scale_for(eig, c: float) -> float:
    """Shared bulk-median noise-scale estimate (detection/tuner agreement)."""
    from de_formulas import estimate_noise_scales

    return estimate_noise_scales(eig[0], c)[0]


def est_eb_spectral(Xc, Ys, eig, rng, tau_override: float | None = None,
                    kmax_hint: int | None = None):
    """SEB soft-trim estimator (ours). Returns (beta_hat, info).

    kmax_hint bounds how many sample eigenvalues may be read as supercritical
    spikes (V5 r-misspecification cells consume cfg.r + cfg.r_misspec).
    """
    n, p = Xc.shape
    c = p / n
    d, V = eig
    l_hat = estimate_spike_profile(d, c, kmax=kmax_hint or 10)
    g2_hat = estimate_gamma2(Xc, Ys, eig, l_hat, c)
    se2 = _noise_scale_for(eig, c)
    if tau_override is not None:
        tau = tau_override
    else:
        taus = [f * float(np.median(d)) for f in EB_TAU_FRAC_GRID]
        vals = [
            seb_predicted_mse(l_hat, g2_hat, d, t, c, n, se2) for t in taus
        ]
        tau = taus[int(np.argmin(vals))]
    w = trim_weights(d, tau)
    beta = spectral_fit(d, V, Xc.T @ Ys / n, w)
    info = {
        "tau": tau,
        "l_hat": l_hat.tolist(),
        "g2_hat": g2_hat.tolist(),
        "grid": len(EB_TAU_FRAC_GRID),
    }
    return beta, info


def est_eb_cv_tau(Xc, Ys, eig, rng):
    """Ablation no-EB: same soft-trim family, tau by EXACT LOO prediction-CV
    (loo_scores_spectral with w = min(1, tau/d)). Isolates what the SEB
    causal objective adds over plain predictive tuning."""
    n, p = Xc.shape
    d, V = eig
    U = left_factors(Xc, eig)
    taus = [fr * float(np.median(d)) for fr in EB_TAU_FRAC_GRID]
    scores = loo_scores_spectral(
        U, Ys, d, lambda t: trim_weights(d, t), taus
    )
    tau = float(taus[int(np.argmin(scores))])
    w = trim_weights(d, tau)
    return spectral_fit(d, V, Xc.T @ Ys / n, w), {"tau": tau}


# ---------------------------------------------------------------------------
# oracle variants (diagnostic upper bounds only; plan Section 8.3)
# ---------------------------------------------------------------------------


def make_oracle_estimator(Q: np.ndarray, gam: np.ndarray):
    """oracle_gamma: beta_hat = beta_ols - (Q Q' beta_ols + Q gamma), i.e.,
    OLS minus the true loading-subspace bias contribution. Upper bound for
    what any spectral method could achieve on the confounding component.
    """

    def f(Xc, Ys, eig, rng=None):
        b_ols = fit_min_norm(Xc, Ys, eig)
        return b_ols - Q @ (Q.T @ b_ols) - Q @ gam

    return f


def eb_oracle_tau_factory(l_true: np.ndarray, gam_true: np.ndarray):
    """eb_oracle_tau: same weight family and objective as est_eb_spectral,
    but fed the TRUE (l, gamma^2) instead of plug-in estimates.

    Frozen definition (preregistration): the oracle ablation isolates tuner
    ESTIMATION error, not objective mismatch. The objective is
    seb_predicted_mse with true parameters; variance/noise scale still uses
    the observable bulk estimate (_noise_scale_for) so only the confounding
    geometry is oracular.
    """

    def f(Xc, Ys, eig, rng=None):
        n, p = Xc.shape
        c = p / n
        d, V = eig
        se2 = _noise_scale_for(eig, c)
        g2_true = np.asarray(gam_true, float) ** 2
        taus = [fr * float(np.median(d)) for fr in EB_TAU_FRAC_GRID]
        vals = [seb_predicted_mse(l_true, g2_true, d, t, c, n, se2) for t in taus]
        tau = taus[int(np.argmin(vals))]
        w = trim_weights(d, tau)
        return spectral_fit(d, V, Xc.T @ Ys / n, w)

    return f
