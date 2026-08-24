"""SCF Phase 1 simulator: DGP generation, estimators, spectral statistics.

Schema: one parquet row per config x rep x estimator tag (research plan
Section 10.3). Estimator tags used by the pilot:
  ols, ridge_fixed (one row per lambda), pca_oracle (k=r), pca_onatski,
  and, for cells with twin_gamma0=True, the gamma=0 twins ols_g0/ridge_fixed_g0
  fitted on the SAME (Q, beta, f, u, eps) draw.

All randomness flows from np.random.SeedSequence([GLOBAL_SEED, cfg_hash, rep])
so runs are reproducible and cells are resumable. Single-thread BLAS must be
set by the launcher before numpy import (research plan Section 10.1).
"""
from __future__ import annotations

import hashlib
import math
import time
from dataclasses import dataclass

import numpy as np

from de_formulas import (
    bbp_location,
    bgn_overlap,
    onatski_select,
    tw_mu_sigma,
    tw_threshold,
)

GLOBAL_SEED = 20260823


@dataclass(frozen=True)
class Config:
    n: int
    p: int
    r: int
    l: tuple[float, ...]
    theta: float
    g: float = 1.0
    sigma_u: float = 1.0
    sigma_eps: float = 1.0
    beta_kind: str = "dense"     # "dense" (A4a) | "aligned" (rung 4)
    conf_kind: str = "dense"     # "dense" | "sparse" (rung 4 sparse loadings)
    loading_kind: str = "gauss"  # "gauss" | "rademacher_half" (WP 2.4 V2)
    error_law: str = "gaussian"  # "gaussian" | "t5" (WP 2.4 robustness)
    hetero_u: bool = False       # row-heteroskedastic u (WP 2.4)
    corr_factors: bool = False   # correlated f (WP 2.4; breaks Var(f)=I)
    r_misspec: int = 0           # r perturbation consumed by priors (V5)
    m2_treatment: bool = False   # generate D block (WP 2.4 / Phase 3 prep)
    m2_tau: float = 1.0          # true scalar treatment coefficient
    delta_g: float = 0.3         # ||delta|| for M2 weak-treatment alignment
    twin_gamma0: bool = False    # run gamma=0 arm on common seeds
    q_fixed: bool = False        # draw loading directions once per config
    reps: int = 200
    profile: str = ""
    label: str = ""

    @property
    def c(self) -> float:
        return self.p / self.n

    @property
    def cid(self) -> str:
        payload = repr(
            (self.n, self.p, self.r, self.l, round(self.theta, 6), self.g,
             self.sigma_u, self.sigma_eps, self.beta_kind, self.conf_kind,
             self.loading_kind, self.error_law, self.hetero_u,
             self.corr_factors, self.r_misspec,
             self.m2_treatment, round(self.m2_tau, 6), round(self.delta_g, 6),
             self.profile, self.label)
        )
        return hashlib.sha256(payload.encode()).hexdigest()[:12]


def gamma_vector(cfg: Config) -> np.ndarray:
    """gamma = g * dir(theta) in factor coordinates (model card Section 2)."""
    gam = np.zeros(cfg.r)
    gam[0] = np.cos(cfg.theta)
    if cfg.r > 1:
        gam[1] = np.sin(cfg.theta)
    return cfg.g * gam


def _rng(cfg: Config, rep: int) -> np.random.Generator:
    return np.random.default_rng(
        np.random.SeedSequence(entropy=[GLOBAL_SEED, int(cfg.cid[:8], 16), rep])
    )


def gen_data(cfg: Config, rep: int):
    """Draw Q, beta, f, u, eps; return dict with X, Y and ingredients.

    Phase 2 extensions (preregistration): sparse-confounding loadings
    (conf_kind="sparse": each loading column supported on ceil(sqrt(p))
    coordinates, rescaled so Lambda'Lambda = diag(sigma_u^2 l_j) exactly),
    aligned beta (beta_kind="aligned": cosine 0.9 with the top factor
    direction), t5 errors, row-heteroskedastic u, correlated factors, and an
    M2 treatment block D = pi'X + delta'f + nu.
    """
    rng = _rng(cfg, rep)
    p, n, r = cfg.p, cfg.n, cfg.r
    if cfg.q_fixed:
        qrng = np.random.default_rng(
            np.random.SeedSequence(entropy=[GLOBAL_SEED, int(cfg.cid[:8], 16), 2 ** 32])
        )
        A = qrng.standard_normal((p, r))
    else:
        A = rng.standard_normal((p, r))
    Q, _ = np.linalg.qr(A)
    svals = cfg.sigma_u * np.sqrt(np.asarray(cfg.l))
    if cfg.conf_kind == "sparse":
        # DISJOINT supports (r * ceil(sqrt(p)) <= p in all grids) so that
        # Lambda'Lambda = diag(sigma_u^2 l_j) holds exactly and Ubasis stays
        # orthonormal; rescaled columns preserve the spike strengths.
        ksup = max(int(round(math.sqrt(p))), r)
        Ldraw = qrng if cfg.q_fixed else rng
        perm = Ldraw.permutation(p)
        cols = []
        for j in range(r):
            support = perm[j * ksup:(j + 1) * ksup]
            col = np.zeros(p)
            col[support] = Ldraw.standard_normal(ksup)
            cols.append(col)
        Acol = np.column_stack(cols)
        norms = np.linalg.norm(Acol, axis=0)
        Lam = Acol * (svals / np.maximum(norms, 1e-12))[None, :]
    elif cfg.loading_kind == "rademacher_half":
        # WP 2.4 V2: Bernoulli-Rademacher loadings on a half support,
        # rescaled so Lambda'Lambda = diag(sigma_u^2 l_j) exactly.
        Ldraw = qrng if cfg.q_fixed else rng
        signs = Ldraw.choice([-1.0, 1.0], size=(p, r))
        support_mask = Ldraw.random((p, r)) < 0.5
        Acol = np.where(support_mask, signs, 0.0)
        norms = np.linalg.norm(Acol, axis=0)
        norms = np.maximum(norms, 0.05 * math.sqrt(p / 2))  # avoid degenerate cols
        norms = np.maximum(norms, 1e-6)
        Lam = Acol * (svals / norms)[None, :]
    else:
        Lam = Q * svals  # p x r with scaled columns

    if cfg.beta_kind == "aligned":
        w = rng.standard_normal(p)
        w -= Q @ (Q.T @ w)
        nw = np.linalg.norm(w)
        w = w / nw if nw > 1e-12 else w
        beta = 0.9 * Q[:, 0] + math.sqrt(1.0 - 0.81) * w
    else:
        beta = rng.standard_normal(p)
        beta /= np.linalg.norm(beta)

    if cfg.corr_factors:
        OmV, _ = np.linalg.qr(rng.standard_normal((r, r)))
        om_evals = rng.uniform(0.5, 1.5, size=r)
        Omega = (OmV * om_evals) @ OmV.T
        f = rng.multivariate_normal(np.zeros(r), Omega, size=n)
    else:
        f = rng.standard_normal((n, r))

    U_raw = rng.standard_normal((n, p))
    if cfg.error_law == "t5":
        # unit-variance t: t_5 / sqrt(5/3)
        U_raw = rng.standard_t(5, size=(n, p)) / math.sqrt(5.0 / 3.0)
        eps_core = rng.standard_t(5, size=n) / math.sqrt(5.0 / 3.0)
    else:
        eps_core = rng.standard_normal(n)
    if cfg.hetero_u:
        row_scale = np.sqrt((1.0 + rng.chisquare(1, size=(n, 1))) / 2.0)
        U_raw = U_raw * row_scale
    U = cfg.sigma_u * U_raw
    X = f @ Lam.T + U
    gam = gamma_vector(cfg)
    noise_eps = cfg.sigma_eps * eps_core
    out = dict(X=X, beta=beta, Q=Q, Lam=Lam, gam=gam, eps=noise_eps)
    if cfg.m2_treatment:
        pi = np.zeros(p)
        kpi = max(3, p // 100)
        pi[:kpi] = 1.0 / math.sqrt(kpi)
        delta = rng.standard_normal(r)
        delta *= cfg.delta_g / max(float(np.linalg.norm(delta)), 1e-12)
        nu_ = cfg.sigma_eps * rng.standard_normal(n)
        D = X @ pi + f @ delta + nu_
        out.update(D=D, pi=pi, delta=delta)
        out["Y"] = cfg.m2_tau * D + X @ beta + f @ gam + noise_eps
    else:
        out["Y"] = X @ beta + f @ gam + noise_eps
    return out


def spectrum(X: np.ndarray):
    """Descending eigenpairs of the sample covariance on the R^p side.

    For p > n the n-side Gram is diagonalized and eigenvectors are mapped back
    via V_j = X' q_j / sqrt(n d_j) (exact normalization).
    """
    n, p = X.shape
    if p <= n:
        d, W = np.linalg.eigh(X.T @ X / n)
    else:
        d, W = np.linalg.eigh(X @ X.T / n)
        d = np.maximum(d, 1e-12)
        W = X.T @ W / np.sqrt(n * d)
    idx = np.argsort(d)[::-1]
    return d[idx], W[:, idx]


def fit_ols(X, Y, eig) -> np.ndarray:
    """OLS (min-norm when p > n): V diag(1/d) V' X'Y/n or its wide analogue."""
    d, V = eig
    n = X.shape[0]
    return V @ ((V.T @ (X.T @ Y / n)) / d)


def fit_ridge(X, Y, eig, lam: float) -> np.ndarray:
    """(Sigmahat + lam I)^{-1} X'Y/n; exact for any p via the shared spectrum."""
    d, V = eig
    n = X.shape[0]
    rhs = V.T @ (X.T @ Y / n)
    return V @ (rhs / (d + lam))


def fit_pca(X, Y, eig, k: int) -> np.ndarray:
    """Regress Y on the top-k sample PCs (k >= 1); k <= 0 gives zero vector."""
    if k <= 0:
        return np.zeros(X.shape[1])
    Vk = eig[1][:, :k]
    delta, *_ = np.linalg.lstsq(X @ Vk, Y, rcond=None)
    return Vk @ delta


def run_rep(cfg: Config, rep: int, lam_grid: tuple[float, ...]):
    """Run one replication; returns (rows, mean_diffs).

    rows: list of dicts (one per estimator tag) with scalar metrics and shared
    spectral statistics.
    mean_diffs: tag -> (beta_hat - beta) vector, accumulated across reps by the
    runner to estimate E[beta_hat] - beta without storing every draw.
    """
    t0 = time.perf_counter()
    data = gen_data(cfg, rep)
    X, Y, beta, Q = data["X"], data["Y"], data["beta"], data["Q"]
    n, p = X.shape
    eig = spectrum(X)
    d, V = eig
    mu_np, sig_np = tw_mu_sigma(n, p)
    lam_max = float(d[0])

    stats_common = {
        "lam_max_cov": lam_max,
        "tw_stat": float((lam_max * n - mu_np) / sig_np),
        "outlier99": bool(lam_max > tw_threshold(n, p, cfg.sigma_u)),
        "bbp_pred": bbp_location(max(cfg.l), cfg.c, cfg.sigma_u),
        "xi1_pred": bgn_overlap(max(cfg.l), cfg.c),
    }
    overlaps = []
    for j in range(min(cfg.r, V.shape[1])):
        overlaps.append(float((V[:, j] @ Q[:, j]) ** 2))

    arms = [("ols", np.nan)]
    for lam in lam_grid:
        arms.append((f"ridge_fixed|{lam}", lam))

    def make_rows(suffix: str, Y_arm: np.ndarray):
        out = []
        diffs = {}
        bh_ols = fit_ols(X, Y_arm, eig)
        ests = [(f"ols{suffix}", np.nan, bh_ols)]
        for _, lam in arms[1:]:
            ests.append((f"ridge_fixed{suffix}", lam, fit_ridge(X, Y_arm, eig, lam)))
        k_on = max(onatski_select(d), 0)
        ests.append((f"pca_onatski{suffix}", np.nan, fit_pca(X, Y_arm, eig, k_on)))
        for tag, lam, bh in ests:
            diff = bh - beta
            row = {
                "config_id": cfg.cid,
                "rep": rep,
                "seed": GLOBAL_SEED,
                "estimator": tag,
                "lambda": lam,
                "k_select": k_on if tag.startswith("pca") else -1,
                "rel_err": float(np.linalg.norm(diff)),
                "runtime_s": 0.0,
            }
            row.update(stats_common)
            for j in range(cfg.r):
                ov = overlaps[j] if j < len(overlaps) else np.nan
                row[f"overlap{j + 1}"] = ov
            out.append(row)
            key = tag if np.isnan(lam) else f"{tag}@{lam}"
            diffs[key] = diff.astype(np.float64)
        return out, diffs

    rows, mean_diffs = make_rows("", Y)

    if cfg.twin_gamma0:
        Y0 = X @ beta + data["eps"]  # same seeds: identical everything but gamma link
        rows0, diffs0 = make_rows("_g0", Y0)
        rows.extend(rows0)
        mean_diffs.update(diffs0)

    runtime = time.perf_counter() - t0
    for r_ in rows:
        r_["runtime_s"] = runtime / len(rows)
    return rows, mean_diffs
