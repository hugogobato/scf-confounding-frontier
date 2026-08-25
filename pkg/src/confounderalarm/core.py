"""confounderalarm: calibrated alarm for hidden dense confounding.

Public API:
    fit_alarm(Y, X, D=None, ...) -> AlarmReport (dict subclass)

The alarm tests, relative to the A4a near-orthogonality ledger of the SCF
model card, whether the response carries a dense-factor link that the
leading design coordinates can see. Outputs: verdict with permutation
p-value, estimated (r_hat, l_hat, c), placement relative to the predicted
detectability frontier g*, a blind-region certificate, and, when a treatment
block D is supplied, the recommended Onatski hard-trim adjustment.
"""
from __future__ import annotations

import numpy as np

from ._spectra import (
    center_columns,
    center_vec,
    noise_floor_bench,
    onatski_select,
    onatski_trim_tau,
    predicted_g_star,
    raw_z_coords,
    spectrum,
    tw_threshold,
    ucm_rho_proxy,
    js_asymmetry,
)

__all__ = ["fit_alarm", "AlarmReport"]
__version__ = "0.1.0"


class AlarmReport(dict):
    """dict with attribute access for convenience."""

    def __getattr__(self, k):
        try:
            return self[k]
        except KeyError as e:
            raise AttributeError(k) from e


DEFAULT_SEED = 20260823


def fit_alarm(Y, X, D=None, alpha: float = 0.05, n_perm: int = 400,
              seed: int = DEFAULT_SEED, kmax: int = 10) -> AlarmReport:
    """Fit the confounding alarm on one dataset.

    Parameters
    ----------
    Y : array (n,), response.
    X : array (n, p), covariate design.
    D : optional array (n,), observed treatment. When given, the report adds
        the trim-then-regress adjustment (Onatski hard trim).
    alpha : test level for the verdict.
    n_perm : permutations used for the null scales and the p-value.
    seed : RNG seed (permutation stream).

    Returns
    -------
    AlarmReport with keys:
        alarm, p_value, statistic T, threshold (permutation q95),
        n, p, c, r_hat, l_hat (top-k, normalized), se2, ktop,
        tw_stat_lammax, outlier99_white,
        g_star (predicted frontier; inf = no supercritical-aligned mass),
        placement in {"below", "near", "above", "blind"},
        certificate (human-readable blind-region statement),
        ucm_rho, js_asym (diagnostics, uncalibrated),
        adjustment (None or dict) when D is given.
    """
    Y = np.asarray(Y, float).ravel()
    X = np.asarray(X, float)
    n, p = X.shape
    Xc = center_columns(X)
    Yc = center_vec(Y)
    eig = spectrum(Xc)
    d = eig[0]
    c = p / n
    se2 = noise_floor_bench(d)
    l_hat = np.maximum(d[:10] / max(se2, 1e-12) - 1.0, 0.0)
    r_hat = int(onatski_select(d))
    ktop = int(min(max(r_hat, 1), kmax))

    mu_np, sig_np = tw_mu_sigma_safe(n, p)
    lam_max = float(d[0])
    tw_stat = float((lam_max * n - mu_np) / sig_np)
    outlier99 = bool(lam_max > tw_threshold(n, p, 1.0))

    # ---- permutation calibration (marginals preserved) -------------------
    rng = np.random.default_rng(seed)
    zeta_obs = raw_z_coords(Xc, Yc, eig, ktop)
    Zperm = np.empty((int(n_perm), ktop))
    rho_perm = np.empty(int(n_perm))
    js_perm = np.empty(int(n_perm))
    for i in range(int(n_perm)):
        yp = rng.permutation(Yc)
        Zperm[i] = raw_z_coords(Xc, yp, eig, ktop)
        rho_perm[i] = ucm_rho_proxy(Xc, yp, eig, l_hat[:ktop], c)
        js_perm[i] = js_asymmetry(Xc, yp, eig, K=ktop + 2)
    s = np.sqrt((Zperm ** 2).mean(axis=0))
    T_perm = np.max(np.abs(Zperm) / np.maximum(s, 1e-12)[None, :], axis=1)
    T_obs = float(np.max(np.abs(zeta_obs) /
                         np.maximum(s, 1e-12)))
    mc95 = float(np.quantile(T_perm, 0.95))
    p_value = float((1 + int((T_perm >= T_obs).sum())) / (n_perm + 1))

    # ---- frontier placement ---------------------------------------------
    dirv = np.ones(ktop) / np.sqrt(ktop)
    g_star = predicted_g_star(l_hat[:ktop], s, dirv, c, n, se2,
                              mc95=mc95, seed=seed + 1)
    if not np.isfinite(g_star):
        placement, cert = "blind", (
            "no supercritical-aligned link is visible to this statistic "
            "class at any g: non-detection here is UNINFORMATIVE about "
            "dense confounding below the design's visibility boundary")
    elif T_obs >= mc95 and g_star <= 1.0 * _g_scale_hint(T_obs, s, l_hat,
                                                         dirv, c, n, se2,
                                                         mc95, seed):
        placement, cert = "above", (
            "alarm fired and the data sit above the predicted detection "
            "frontier; non-detection would have been informative, so treat "
            "the alarm as evidence of a dense factor link")
    elif p_value <= alpha:
        placement, cert = "above", (
            "alarm fired at level alpha")
    elif T_obs > 0.5 * mc95:
        placement, cert = "near", (
            "statistic sits within a factor ~2 of the calibrated threshold; "
            "treat the verdict as marginal")
    else:
        placement, cert = "below", (
            "below the calibrated threshold AND below the predicted "
            "frontier: non-detection is uninformative about harmful dense "
            "confounding in the blind region (bias can be O(1) while no "
            "second-moment alarm fires)")

    rho_obs = float(ucm_rho_proxy(Xc, Yc, eig, l_hat[:ktop], c))
    js_obs = float(js_asymmetry(Xc, Yc, eig, K=ktop + 2))

    report = AlarmReport(
        alarm=bool(p_value <= alpha),
        p_value=p_value,
        statistic=T_obs,
        threshold=mc95,
        n=n, p=p, c=float(c),
        r_hat=r_hat, ktop=ktop,
        l_hat=[float(v) for v in l_hat[:ktop]],
        se2=float(se2),
        tw_stat_lammax=tw_stat,
        outlier99_white=outlier99,
        g_star=(None if not np.isfinite(g_star) else float(g_star)),
        placement=placement,
        certificate=cert,
        ucm_rho=rho_obs,
        js_asym=js_obs,
        ucm_p=float((1 + int((rho_perm >= rho_obs).sum())) / (n_perm + 1)),
        js_p=float((1 + int((js_perm >= js_obs).sum())) / (n_perm + 1)),
        adjustment=None,
        alpha=alpha, n_perm=int(n_perm), seed=int(seed),
    )
    if D is not None:
        Dc = center_vec(np.asarray(D, float).ravel())
        tau_trim, tau_ols, k_trim = onatski_trim_tau(Dc, Xc, Yc, eig,
                                                     r_hint=r_hat)
        report["adjustment"] = {
            "method": "onatski_trim",
            "k": int(k_trim),
            "tau_trim": float(tau_trim),
            "tau_ols": float(tau_ols),
            "note": ("hard-trim adjustment per the Phase-2 M2 result; "
                     "tuned soft weights are NOT recommended"),
        }
    return report


def _g_scale_hint(T_obs, s, l_hat, dirv, c, n, se2, mc95, seed) -> float:
    """Crude data-driven g estimate used only to label 'above' placement."""
    sup = [j for j in range(len(dirv)) if dirv[j] != 0]
    num = 0.0
    den = 0.0
    for j in sup:
        from ._spectra import bbp_location

        dj = bbp_location(float(l_hat[j]), c, se2)
        omega = float(min(1.0, (1.0 + l_hat[j]) / (c + l_hat[j]))) if c > 1 \
            else 1.0
        sj = max(float(s[j]), 1e-12)
        slope = (np.sqrt(n) * omega * np.sqrt(se2 * l_hat[j]) * abs(dirv[j]) /
                 (np.sqrt(dj) * sj))
        num += slope * (T_obs / max(mc95, 1e-9))
        den += slope ** 2
    return float(num / max(den, 1e-12))


def tw_mu_sigma_safe(n, p):
    from ._spectra import tw_mu_sigma

    return tw_mu_sigma(n, p)
