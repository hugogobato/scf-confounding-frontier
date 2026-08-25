"""T1 capture-law theorem tests (Route A-prime, closed 2026-08-25).

Falsifiers for docs/theory_T1_capture_law.md:
  - corrected assembly of the two exact channel pieces ch_a/ch_b
  - collapse law sum -> sqrt(l)/(c+l)
  - artifact theorem R2: E[q'Pi q] -> (c-1)/(c+l)
  - conditional directional mean incl. the fit-artifact channel
  - r=2 decoupling
  - T1.c ridge capture (proved form) vs simulation + collapse checks

All statistics are computed with O(n^2) probe algebra (no p x p projector).
Single-threaded BLAS; deterministic seeds.
"""

import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "code"))

from de_formulas import ridge_capture, ridge_capture_superseded  # noqa: E402


def channel_pieces(l, c, n_reps, n, seed):
    """Per-rep exact pieces ch_a, ch_b, R2 for r=1 min-norm OLS."""
    rng = np.random.default_rng(seed)
    p = int(round(c * n))
    out = np.empty((n_reps, 3))
    for i in range(n_reps):
        f = rng.standard_normal(n)
        U = rng.standard_normal((n, p))
        q = rng.standard_normal(p)
        q /= np.linalg.norm(q)
        eps = rng.standard_normal(n) * 0.7
        u = U @ q
        G = (l * np.outer(f, f) + np.sqrt(l) * (np.outer(f, u) + np.outer(u, f))
             + U @ U.T)
        Ginv = np.linalg.inv(G)
        e = f / np.linalg.norm(f)
        ch_a = np.sqrt(l) * (f @ f) * (e @ Ginv @ e)
        ch_b = np.sqrt((f @ f)) * (u @ Ginv @ e)
        Xq = np.sqrt(l) * f + u
        w = Ginv @ Xq
        XtGXq = np.sqrt(l) * q * (f @ w) + U.T @ w   # X' G^-1 Xq, O(n p + n^2)
        out[i] = (ch_a, ch_b, 1.0 - q @ XtGXq)
    return out


def directional_mean(l_vec, c, gammas, n_reps, n, seed):
    """Fixed (Q, beta), fresh (f, U, eps): per-rep <beta_hat-beta, q_j>."""
    rng = np.random.default_rng(seed)
    p = int(round(c * n))
    r = len(l_vec)
    Q, _ = np.linalg.qr(np.random.default_rng(seed + 1).standard_normal((p, r)))
    beta = np.random.default_rng(seed + 2).standard_normal(p)
    beta /= np.linalg.norm(beta)
    bQ = Q.T @ beta
    L = np.sqrt(np.asarray(l_vec))
    acc = np.zeros(r)
    for _ in range(n_reps):
        F = rng.standard_normal((n, r))
        U = rng.standard_normal((n, p))
        eps = rng.standard_normal(n) * 0.7
        X = (F * L[None, :]) @ Q.T + U
        Y = X @ beta + F @ gammas + eps
        G = X @ X.T
        Ginv = np.linalg.inv(G)
        d = X.T @ (Ginv @ Y) - beta
        for j in range(r):
            acc[j] += d @ Q[:, j]
    return acc / n_reps, bQ


TOL_REL_PCT = 5.0
TOL_ABS_R2 = 0.03
TOL_ABS_DIR = 0.04


def test_channel_assembly_and_collapse():
    """ch_a, ch_b piece formulas and the collapse sqrt(l)/(c+l), t != 1 cells."""
    cells = [(4.0, 2.0), (6.708, 5.0), (16.0, 8.0)]
    for l, c in cells:
        res = channel_pieces(l, c, n_reps=200, n=300, seed=int(1000 * l + 10 * c))
        t = 1.0 / (c - 1.0)
        pa = np.sqrt(l) * t * (1 + t) / (1 + t * (1 + l))
        pb = -np.sqrt(l) * t * t / (1 + t * (1 + l))
        tot = res[:, :2].mean(axis=0).sum()
        target = np.sqrt(l) / (c + l)
        dev_pct = abs(tot - target) / abs(target) * 100
        assert dev_pct < TOL_REL_PCT, f"sum dev {dev_pct:.2f}% at ({l},{c})"
        # individual pieces at tighter tolerance in absolute terms
        assert abs(res[:, 0].mean() - pa) < 0.05, f"ch_a at ({l},{c})"
        assert abs(res[:, 1].mean() - pb) < 0.05, f"ch_b at ({l},{c})"


def test_r2_artifact_theorem():
    for l, c in [(4.0, 2.0), (6.708, 5.0), (2.0, 1.5)]:
        res = channel_pieces(l, c, n_reps=200, n=300, seed=int(2000 * l + 10 * c))
        pred = (c - 1.0) / (c + l)
        assert abs(res[:, 2].mean() - pred) < TOL_ABS_R2, f"R2 at ({l},{c})"


def test_conditional_directional_mean_with_artifact():
    l, c, g = 6.708, 5.0, 1.3
    dm, bq = directional_mean([l], c, [g], n_reps=250, n=300, seed=42)
    cap = (1 + l) / (c + l)
    pred = (cap - 1) * bq[0] + cap * np.sqrt(l) / (1 + l) * g
    assert abs(dm[0] - pred) < TOL_ABS_DIR


def test_r2_decoupling():
    l_vec, c = np.array([6.708, 0.8]), 5.0
    g = np.array([1.3, 0.9])
    dm, bq = directional_mean(l_vec, c, g, n_reps=250, n=300, seed=777)
    for j in range(2):
        cap = (1 + l_vec[j]) / (c + l_vec[j])
        pred = ((cap - 1) * bq[j]
                + cap * np.sqrt(l_vec[j]) / (1 + l_vec[j]) * g[j])
        assert abs(dm[j] - pred) < TOL_ABS_DIR, f"spike {j}"


def sim_dir_ridge(l, c, gam, lam, n_reps, n, seed):
    """Per-rep <beta_hat_ridge - beta, q> with fixed (Q, beta)."""
    rng = np.random.default_rng(seed)
    p = int(round(c * n))
    Q, _ = np.linalg.qr(np.random.default_rng(seed + 1).standard_normal((p, 1)))
    q = Q[:, 0]
    beta = np.random.default_rng(seed + 2).standard_normal(p)
    beta /= np.linalg.norm(beta)
    bq = float(beta @ q)
    acc = 0.0
    for _ in range(n_reps):
        f = rng.standard_normal(n)
        U = rng.standard_normal((n, p))
        eps = rng.standard_normal(n) * 0.7
        X = np.sqrt(l) * np.outer(f, q) + U
        Y = X @ beta + gam * f + eps
        S = X.T @ X / n
        br = np.linalg.solve(S + lam * np.eye(p), X.T @ Y / n)
        acc += (br - beta) @ q
    return acc / n_reps, bq


def test_ridge_capture_corrected():
    """T1.c: proved shifted-resolvent ridge capture vs simulation."""
    g = 1.3
    for l, c in [(4.0, 2.0), (0.5, 2.0)]:
        for lam in [0.25, 2.0]:
            dm, bq = sim_dir_ridge(l, c, g, lam, n_reps=150, n=300,
                                   seed=int(100 * l + 7 * lam))
            cap = float(ridge_capture(np.array([l]), lam, c)[0])
            pred = (cap - 1.0) * bq + cap * np.sqrt(l) / (1.0 + l) * g
            assert abs(dm - pred) < TOL_ABS_DIR, f"ridge at ({l},{c},{lam})"


def test_ridge_collapse_checks():
    """lam->0 recovers min-norm capture; lam->inf -> 0."""
    l, c = 6.708, 5.0
    cap0 = float(ridge_capture(np.array([l]), 1e-10, c)[0])
    assert abs(cap0 - (1.0 + l) / (c + l)) < 1e-3
    cap_inf = float(ridge_capture(np.array([l]), 1e8, c)[0])
    assert cap_inf < 1e-3


if __name__ == "__main__":
    import sys

    sys.exit(pytest.main([__file__, "-v"]))
