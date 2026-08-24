"""SCF Phase 2 detection statistics (WP 2.3). Frozen spec:
docs/detection_statistics.md. All statistics consume the shared spectrum
(descending eigenpairs of X_c'X_c/n) plus the standardized response.

Per-rep outputs are plain floats collected into parquet rows by the runner;
thresholds follow the frozen spec (analytic variants computed here,
MC-calibrated variants computed at analysis time from pooled null reps).
"""
from __future__ import annotations

import numpy as np

from de_formulas import (
    TW1_Q95,
    TW1_Q99,
    maxz_null_var as df_maxz_null_var,
    onatski_select,
    tw_threshold,
    tw_width_cov_scale,
)


def ktop_default(eig, cap: int = 10) -> int:
    """Number of spike coordinates used by S2/B1: Onatski selection, min 1."""
    return int(min(max(onatski_select(eig[0]), 1), cap))


def compute_stats(
    Xc: np.ndarray,
    Yraw: np.ndarray,
    eig,
    sigma_u: float = 1.0,
    ktop: int | None = None,
):
    """Yraw is the centered raw response; it is standardized internally
    (scale-invariant statistics)."""
    sd_y = float(np.std(Yraw)) or 1.0
    Ys = Yraw / sd_y
    """Return dict of detection statistics for one dataset.

    Keys:
      lam_max_cov, tw_stat, outlier99          -> S0 scree baseline
      t_aug, s1_thresh_analytic                -> S1 aug_bbp
      t_maxz, s2_thresh_bonf                   -> S2 maxz_cal
      lss_rank1                                -> S3 reduction (rank-one LSS)
      f_pcs                                    -> B1 partial F on PC scores
      b_norm2, z_coords                        -> probe features / diagnostics
      ktop                                     -> coordinate count used
    """
    n, p = Xc.shape
    c = p / n
    d, V = eig
    b = Xc.T @ Ys / n
    if ktop is None:
        ktop = ktop_default(eig)

    # ---- S0: scree / TW99 visibility -------------------------------------
    mu_np = None
    lam_max = float(d[0])
    from de_formulas import tw_mu_sigma

    mu_np, sig_np = tw_mu_sigma(n, p)
    tw_stat = float((lam_max * n - mu_np) / sig_np)
    outlier99 = bool(lam_max > tw_threshold(n, p, sigma_u))

    # ---- S2: calibrated max-z over spike coordinates ----------------------
    # F12 (Erratum 1): Var(z_j) = (sigma_eps2 + d_j/c)/sigma_y2 with
    # sigma_y2 = tr(Sigmahat)/p + sigma_eps2; scales from shared helper.
    kk = min(ktop, len(d))
    w = V[:, :kk].T @ b
    z = np.sqrt(n) * w / np.sqrt(d[:kk])
    from de_formulas import estimate_noise_scales

    se2_hat, sigma_y2 = estimate_noise_scales(d, c)
    var_cal = np.array([
        df_maxz_null_var(dj, c, se2_hat, sigma_y2) for dj in d[:kk]
    ])
    t_maxz = float(np.max(np.abs(z) / np.sqrt(var_cal)))
    from scipy.stats import norm

    s2_thresh = float(norm.ppf(1.0 - 0.05 / (2 * kk)))

    # ---- S1: augmented secular root ---------------------------------------
    w_all = V.T @ b
    dpos = np.maximum(d, 1e-300)
    lo = float(d[0]) * (1.0 + 1e-9) + 1e-15
    s_total = float(np.sum(w_all ** 2))
    disc = (1.0 + lo) ** 2 - 4.0 * (lo + s_total)
    hi = 0.5 * ((1.0 + lo) + np.sqrt(max(disc, 1e-300))) * (1.0 + 1e-6) + 1e-6

    def f_aug(lam):
        return lam - 1.0 - float(np.sum(w_all ** 2 / (lam - dpos)))

    flo = f_aug(lo)
    if flo > 0:  # degenerate: root pinned at top edge
        t_aug = lo
    else:
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if f_aug(mid) > 0:
                hi = mid
            else:
                lo = mid
            if hi - lo < 1e-11 * max(1.0, hi):
                break
        t_aug = 0.5 * (lo + hi)

    # analytic H0 threshold: plug-in null root + TW95 width at dims (n, p+1)
    from de_formulas import aug_null_root

    lam0 = aug_null_root(d, max(sigma_y2, 1e-6))
    width = tw_width_cov_scale(n, p + 1, sigma_u) * max(sigma_y2, 1.0) ** 0.5
    s1_thresh = float(lam0 + TW1_Q95 * width)

    # ---- S3 reduction: rank-one LSS ----------------------------------------
    lss_rank1 = float(n * np.sum(w_all ** 2 / dpos))

    # ---- B1: partial F of Y on top-ktop PC scores ---------------------------
    S = Xc @ V[:, :kk]
    rss_r = float(np.sum(Ys ** 2))
    coef, *_ = np.linalg.lstsq(S, Ys, rcond=None)
    rss_f = float(np.sum((Ys - S @ coef) ** 2))
    denom = rss_f / max(n - kk, 1)
    f_pcs = ((rss_r - rss_f) / kk) / max(denom, 1e-300)

    return {
        "lam_max_cov": lam_max,
        "tw_stat": tw_stat,
        "outlier99": outlier99,
        "t_aug": float(t_aug),
        "s1_thresh_analytic": s1_thresh,
        "lam0_plugin": float(lam0),
        "t_maxz": t_maxz,
        "s2_thresh_bonf": s2_thresh,
        "lss_rank1": lss_rank1,
        "f_pcs": float(f_pcs),
        "b_norm2": float(np.sum(b ** 2)),
        "z_top": float(z[0]),
        "ktop": int(kk),
    }


def rejections(stats: dict) -> dict:
    """Frozen decision rules applied to a stats dict."""
    return {
        "rej_s0_tw99": bool(stats["outlier99"]),
        "rej_s1_analytic": bool(stats["t_aug"] > stats["s1_thresh_analytic"]),
        "rej_s2_bonf": bool(stats["t_maxz"] > stats["s2_thresh_bonf"]),
        "rej_b1_f95": bool(stats["f_pcs"] > _f_ppf95(stats["ktop"], 10_000)),
    }


_F95_CACHE: dict[tuple[int, int], float] = {}


def mc_thresholds(null_rows: list[dict], q: float = 0.95) -> dict[str, float]:
    """Monte-Carlo decision thresholds from pooled null reps (frozen rule).

    The gate statistics are the MC-calibrated rejections; the analytic
    variants are co-recorded but reported as-is (v1 analytic widths are
    known-miscalibrated for the augmented deformation: measured sd/width
    ratio ~ 8.5 at n=600 before any sweep data, see Erratum notes in
    docs/detection_statistics.md).
    """
    t_aug = np.array([r["t_aug"] for r in null_rows])
    t_z = np.array([r["t_maxz"] for r in null_rows])
    lss = np.array([r["lss_rank1"] for r in null_rows])
    f_pcs = np.array([r["f_pcs"] for r in null_rows])
    return {
        "q_t_aug": float(np.quantile(t_aug, q)),
        "q_t_maxz": float(np.quantile(t_z, q)),
        "q_lss": float(np.quantile(lss, q)),
        "q_f": float(np.quantile(f_pcs, q)),
    }


def _f_ppf95(df1: int, df2: int) -> float:
    key = (df1, min(df2, 100000))
    if key not in _F95_CACHE:
        from scipy.stats import f as f_dist

        _F95_CACHE[key] = float(f_dist.ppf(0.95, df1, df2))
    return _F95_CACHE[key]


# ---------------------------------------------------------------------------
# Numerical Le Cam probe helpers (frozen feature map; lazy sklearn import)
# ---------------------------------------------------------------------------


def probe_features(eig, stats: dict, r_hint: int = 4) -> np.ndarray:
    """Frozen feature map (docs/detection_statistics.md): 12 + 4 values."""
    d = eig[0]
    logs = np.log(np.maximum(d[:10], 1e-12))
    if len(logs) < 10:
        logs = np.pad(logs, (0, 10 - len(logs)))
    extra = [
        float(np.sum(d)),
        stats["t_aug"],
        stats["b_norm2"],
        stats["z_top"],
    ]
    feats = np.concatenate([logs, np.array(extra)])
    if r_hint > 4:
        feats = np.concatenate([feats, np.zeros(4)])
    return feats


def lecam_auc(H0_feats: np.ndarray, H1_feats: np.ndarray, seed: int = 0):
    """GBM + median-heuristic MMD two-sample probe; returns (auc_gbm, auc_mmd).

    Computational probe ONLY (not information-theoretic): declaration
    threshold AUC <= 0.55 for both probes per frozen spec. 50/50 held-out
    split, stratified; GBM = HistGradientBoostingClassifier defaults.
    """
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import train_test_split

    X = np.vstack([H0_feats, H1_feats])
    y = np.concatenate([np.zeros(len(H0_feats)), np.ones(len(H1_feats))])
    Xtr, Xte, ytr, yte = train_test_split(
        X, y, test_size=0.5, random_state=seed, stratify=y
    )
    clf = HistGradientBoostingClassifier(random_state=seed)
    clf.fit(Xtr, ytr)
    auc_gbm = float(roc_auc_score(yte, clf.predict_proba(Xte)[:, 1]))

    mu = Xtr.mean(axis=0)
    sd = Xtr.std(axis=0) + 1e-12
    A = (Xtr[ytr == 1] - mu) / sd
    B = (Xtr[ytr == 0] - mu) / sd
    Z = (Xte - mu) / sd

    def kmat(A_, B_):
        aa = np.sum(A_ * A_, axis=1)[:, None]
        bb = np.sum(B_ * B_, axis=1)[None, :]
        return np.maximum(aa + bb - 2.0 * (A_ @ B_.T), 0.0)

    med = float(np.median(kmat(A[:500], B[:500]))) or 1.0
    gamma = 1.0 / med
    scores = np.exp(-gamma * kmat(Z, A)).mean(axis=1) - np.exp(
        -gamma * kmat(Z, B)
    ).mean(axis=1)
    auc_mmd = float(roc_auc_score(yte, scores))
    return auc_gbm, auc_mmd
