"""Self-contained spectral machinery for confounderalarm.

Conventions (SCF model card): sigma_u^2 estimated from the design's lower
spectrum; population eigenvalues tau_j = 1 + l_j on the normalized scale;
c = p/n. Everything here is pure numpy/scipy so the package has no path
dependence on the research repository.
"""
from __future__ import annotations

import numpy as np

ONATSKI_CRIT = np.array(
    [2.19, 2.09, 2.04, 2.01, 1.99, 1.97, 1.96, 1.95, 1.94, 1.93])
TW1_Q95 = 0.9793
TW1_Q99 = 2.0234


def center_columns(X: np.ndarray) -> np.ndarray:
    return X - X.mean(axis=0, keepdims=True)


def center_vec(y: np.ndarray) -> np.ndarray:
    return y - y.mean()


def spectrum(Xc: np.ndarray):
    """Descending eigenpairs of Xc'Xc/n; wide designs via the n-side Gram."""
    n, p = Xc.shape
    if p <= n:
        d, W = np.linalg.eigh(Xc.T @ Xc / n)
    else:
        d, W = np.linalg.eigh(Xc @ Xc.T / n)
        d = np.maximum(d, 1e-12)
        W = Xc.T @ W / np.sqrt(n * d)
    idx = np.argsort(d)[::-1]
    return d[idx], W[:, idx]


def onatski_select(eigs_desc: np.ndarray, kmax: int = 10) -> int:
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


def tw_mu_sigma(n: int, p: int) -> tuple[float, float]:
    sn, sp = np.sqrt(n - 1.0), np.sqrt(float(p))
    mu = (sn + sp) ** 2
    sigma = (sn + sp) * (1.0 / sn + 1.0 / sp) ** (1.0 / 3.0)
    return mu, sigma


def tw_threshold(n: int, p: int, sigma2: float = 1.0,
                 q: float = TW1_Q99) -> float:
    mu, sigma = tw_mu_sigma(n, p)
    return sigma2 * (mu + q * sigma) / n


def bbp_location(l: float, c: float, sigma2: float = 1.0) -> float:
    if l <= np.sqrt(c):
        return sigma2 * (1 + np.sqrt(c)) ** 2
    return sigma2 * (1 + l) * (l + c) / l


def minnorm_capture(l: np.ndarray, c: float) -> np.ndarray:
    l = np.asarray(l, float)
    return (1.0 + l) / (c + l)


def noise_floor_bench(d_desc: np.ndarray) -> float:
    """Frozen benchmark convention: se2 = max(q25, 1e-3 * mean)."""
    return float(max(np.quantile(d_desc, 0.25),
                     1e-3 * float(np.mean(d_desc))))


def raw_z_coords(Xc: np.ndarray, Yc: np.ndarray, eig, ktop: int):
    """zeta_j = sqrt(n) v_j'(Xc'Yc/n)/sqrt(d_j) on the RAW centered response."""
    n = len(Yc)
    d, V = eig
    k = min(max(int(ktop), 1), len(d))
    b = Xc.T @ Yc / n
    return np.sqrt(n) * (V[:, :k].T @ b) / np.sqrt(d[:k])


def ucm_rho_proxy(Xc, Yc, eig, l_hat: np.ndarray, c: float) -> float:
    """Response-aware confounding-variance share (Rendsburg-et-al.-spirit
    proxy; APPROXIMATE transcription of their estimator)."""
    from scipy.stats import norm  # noqa: F401  (placeholder import guard)

    n, p = Xc.shape
    d, V = eig
    k = min(len(l_hat), 10)
    b = Xc.T @ Yc / n
    sd_y = float(np.std(Yc)) or 1.0
    g2 = []
    for j in range(k):
        wj = float(V[:, j] @ b)
        zj2 = n * wj ** 2 / d[j]
        var_cal = (1.0 + d[j] / c) / max(sd_y ** 2, 1e-12)
        excess = max(0.0, zj2 - var_cal)
        g2.append(excess * (1.0 + l_hat[j]) / (n * max(l_hat[j], 1e-6)))
    conf_var = float(np.sum(np.asarray(g2[:k]) *
                            np.asarray(l_hat[:k]) /
                            (1.0 + np.asarray(l_hat[:k]))))
    return float(conf_var / max(float(np.mean(d)), 1e-12))


def js_asymmetry(Xc, Yc, eig, K: int = 4) -> float:
    """Janzing-Schoelkopf-style spectral asymmetry (APPROXIMATE one-response
    adaptation of janzing2018detecting): relative eigenvalue drops of the
    design covariance after removing the response-explained rank-one part."""
    d, V = eig
    n = len(Yc)
    m = Xc.T @ Yc / n
    sy = float(np.mean(Yc ** 2))
    if sy <= 0:
        return 0.0
    w = V.T @ m
    K = int(min(K, len(d) - 1))

    def h(lam):
        return float(np.sum(w ** 2 / (sy * (d - lam))))

    drops = []
    for i in range(K):
        lo, hi = d[i + 1], d[i]
        if hi - lo < 1e-14 * max(hi, 1e-300):
            continue
        eps_hi = 1e-9 * (hi - lo)
        if h(hi - eps_hi) <= 1.0:
            continue
        a_, b_ = lo + eps_hi, hi - eps_hi
        for _ in range(80):
            mid = 0.5 * (a_ + b_)
            if h(mid) > 1.0:
                b_ = mid
            else:
                a_ = mid
            if b_ - a_ < 1e-11 * hi:
                break
        drops.append((d[i] - 0.5 * (a_ + b_)) / d[i])
    return float(max(drops)) if drops else 0.0


def predicted_g_star(l_hat: np.ndarray, coord_scales: np.ndarray,
                     dirv: np.ndarray, c: float, n: int, se2: float,
                     mc95: float, seed: int = 0) -> float:
    """F12-law frontier: smallest g with predicted alarm power >= 0.8.

    Per-unit-g mean shift of coordinate j (raw-zeta units, standardized by
    the permutation-estimated scale s_j):

        m_j(g)/s_j = g * dir_j * omega_j * sqrt(n se2 l_j) /
                     (sqrt(bbp_location(l_j, c, se2)) * s_j),

    omega_j = clipped min-norm capture at c > 1 else 1; the power curve uses
    a seeded Gaussian-max MC construction (documented approximation).
    """
    sup = [j for j in range(min(len(dirv), len(coord_scales)))
           if dirv[j] != 0]
    if not sup:
        return float("inf")
    slope = np.zeros(len(sup))
    for i, j in enumerate(sup):
        lj = float(l_hat[j])
        dj = bbp_location(lj, c, se2)
        omega = float(np.clip(minnorm_capture(np.array([lj]), c)[0], 0, 1)) \
            if c > 1 else 1.0
        slope[i] = (np.sqrt(n) * omega * np.sqrt(se2 * lj) * abs(dirv[j]) /
                    (np.sqrt(dj) * max(float(coord_scales[j]), 1e-12)))
    if np.max(slope) <= 0:
        return float("inf")
    rr = np.random.default_rng(seed)
    z0 = np.abs(rr.normal(size=(20000, len(sup))))
    thr_sim = float(np.quantile(z0.max(axis=1), 0.95))
    scale = mc95 / thr_sim
    for g in np.linspace(0.01, 20.0, 400):
        z1 = np.abs(rr.normal(size=(20000, len(sup))) +
                    slope[None, :] * g).max(axis=1)
        if float((z1 * scale > mc95).mean()) >= 0.8:
            return round(float(g), 3)
    return float("inf")


def onatski_trim_tau(Dc: np.ndarray, Xc: np.ndarray, Yc: np.ndarray, eig,
                     r_hint: int | None = None):
    """Trim-then-regress treatment coefficient (Frisch-Waugh on top-k PCs).

    k = max(Onatski selection, r_hint, 1); returns (tau_trim, tau_ols, k).
    This is the post-G3 recommended adjustment: hard trim, NOT tuned soft
    weights (Phase 2 estimation kill, pca_onatski dominance 89/97)."""
    d, V = eig
    k = int(max(onatski_select(d), r_hint or 0, 1))
    S = Xc @ V[:, :k]

    def tau(A):
        coef, *_ = np.linalg.lstsq(A, Yc, rcond=None)
        return float(coef[0])

    return tau(np.column_stack([Dc, S])), tau(np.column_stack([Dc, Xc])), k
