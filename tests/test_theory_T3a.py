"""T3(a) augmented-eigenvalue alarm tests (revised closure 2026-08-25,
docs/theory_T3a_eigenvalue_contiguity.md).

The planned OMH-contiguity/impossibility claim was FALSIFIED by these
falsifiers; the assertions below pin the corrected record:

  F-A   frozen blind-strata csv (recorded surrogate-S1 + S2 artifacts)
  id    t_aug-root solver == lambda_max(explicit M_aug)
  BC    population detachment boundary constants (Theorem A, corrected
        normalization sigma_y^2 = A + g^2)
  F-B'  consistent DOWNWARD separation of lambda_max(M_aug) at c=0.2 sub
        (falsifies contiguity)
  F-B"  fixed-point geometry control: upward separation despite flat Q and
        bulk functionals (falsifies the probe-disguise guess)
  F-C   low-edge channel opens at c=5 sub (prediction P-1)
"""

import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pandas as pd
import pytest
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "code"))

ROOT = os.path.join(os.path.dirname(__file__), "..")


def _draw(l, c, n, theta, rng, g):
    """M1 draw, model card conventions, Haar beta per rep."""
    p = int(round(c * n))
    r = len(l)
    Q, _ = np.linalg.qr(rng.standard_normal((p, r)))
    sv = np.sqrt(np.asarray(l))
    beta = rng.standard_normal(p)
    beta /= np.linalg.norm(beta)
    gam = np.zeros(r)
    gam[0] = g * np.cos(theta)
    if r > 1:
        gam[1] = g * np.sin(theta)
    f = rng.standard_normal((n, r))
    U = rng.standard_normal((n, p))
    eps = rng.standard_normal(n)
    X = (f * sv[None, :]) @ Q.T + U
    Y = X @ beta + f @ gam + eps
    return X, Y


def _aug_stats(X, Y):
    """Robust lambda_max / lambda_min of M_aug via secular roots."""
    from simulator import spectrum

    Xc = X - X.mean(0, keepdims=True)
    Yc = Y - Y.mean()
    Ys = Yc / (float(np.std(Yc)) or 1.0)
    d, V = spectrum(Xc)
    b = Xc.T @ Ys / X.shape[0]
    w2 = (V.T @ b) ** 2
    dp = np.maximum(d, 1e-300)

    def f(lam):
        return lam - 1.0 - float(np.sum(w2 / (lam - dp)))

    def root_above(lo):
        hi = lo + float(np.sum(w2)) + 1.0
        for _ in range(120):
            mid = 0.5 * (lo + hi)
            if f(mid) > 0:
                hi = mid
            else:
                lo = mid
        return 0.5 * (lo + hi)

    def root_below(hi):
        lo = hi - float(np.sum(w2)) - 1.0
        for _ in range(120):
            mid = 0.5 * (lo + hi)
            if f(mid) > 0:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    mx = root_above(float(d[0]) * (1 + 1e-9) + 1e-12)
    mn = root_below(float(d[-1]) * (1 - 1e-9) - 1e-12)
    return mx, mn


def _auc_greater(x0, x1):
    u = mannwhitneyu(x1, x0, alternative="greater").statistic
    return float(u / (len(x0) * len(x1)))


def _mx_cell(l, c, n, theta, gs, reps, seed):
    rng = np.random.default_rng(seed)
    out = {g: [] for g in gs}
    mns = {g: [] for g in gs}
    for _ in range(reps):
        for g in gs:
            X, Y = _draw(l, c, n, theta, rng, g)
            mx, mn = _aug_stats(X, Y)
            out[g].append(mx)
            mns[g].append(mn)
    return {g: np.array(v) for g, v in out.items()}, {g: np.array(v) for g, v in mns.items()}


def test_frozen_blind_strata_recorded_artifacts():
    """Recorded artifacts only: in strata with infinite predicted frontier,
    surrogate-S1 and calibrated-S2 power stay <= 0.25 (S0 excluded by
    design - it alarms on design spikes). The S1 column is the DEGENERATE
    SURROGATE (Section 5 erratum); annotated, not re-interpreted."""
    df = pd.read_csv(os.path.join(ROOT, "results", "power_surface.csv"))
    blind = df[df.g_pred_S2.astype(str) == "inf"]
    assert len(blind) >= 27
    for col in ("pow_S2_cal", "pow_S1_cal"):
        assert blind[col].max() <= 0.25, (col, blind[col].max())


def test_secular_root_identity():
    rng = np.random.default_rng(5)
    X, Y = _draw([0.5], 0.5, 120, np.pi / 6, rng, 1.0)
    Xc = X - X.mean(0, keepdims=True)
    Yc = Y - Y.mean()
    Ys = Yc / np.std(Yc)
    from simulator import spectrum

    d, V = spectrum(Xc)
    b = Xc.T @ Ys / X.shape[0]
    M = np.zeros((len(d) + 1, len(d) + 1))
    M[0, 0] = 1.0
    M[0, 1:] = V.T @ b
    M[1:, 0] = M[0, 1:]
    M[1:, 1:] += np.diag(d)
    lam_direct = float(np.linalg.eigvalsh(M)[-1])
    lam_sec, _ = _aug_stats(X, Y)
    assert abs(lam_direct - lam_sec) < 1e-8 * max(1.0, lam_direct)


def test_boundary_constants_theorem_a():
    """Corrected normalization: top never-detach at Phase-2 sub cells;
    bottom wake iff c > 4."""
    ratios = []
    for c in (0.2, 0.8, 2.0):
        sc = np.sqrt(c)
        l = np.full(3, 0.5 * sc)
        dir2 = np.array([np.cos(np.pi / 6) ** 2, np.sin(np.pi / 6) ** 2, 0.0])
        omega_e = float(np.sum(l * dir2 / (c + 2 * sc - l)))
        B = c + 2 * sc
        ratios.append(omega_e / B)
        assert omega_e < B  # never detached at any g
    assert abs(ratios[0] - 0.235) < 0.02 and abs(ratios[2] - 0.036) < 0.005
    for c in (0.2, 0.8, 2.0, 5.0):
        sc = np.sqrt(c)
        l = 0.5 * sc
        omega = float(l * (np.cos(np.pi / 6) ** 2 + np.sin(np.pi / 6) ** 2))
        b_m = (1 - sc) ** 2
        wakes = (1 - b_m) * ((1 + l) - b_m) < omega
        assert wakes == (c > 4.0), (c, wakes)


def test_FB_downward_separation_falsifies_contiguity():
    """F-B': inside the no-detachment region lambda_max(M_aug) separates
    consistently DOWNWARD (AUC << 0.5): spec(M_aug) laws are not contiguous."""
    l = [0.5 * np.sqrt(0.2)] * 3
    res, _ = _mx_cell(l, 0.2, 3200, np.pi / 6, [0.0, 3.2], reps=22, seed=21)
    auc = _auc_greater(res[0.0], res[3.2])
    assert auc < 0.30, auc


def test_FB_fixed_point_geometry_control():
    """F-B": at omega = omega* (theta=pi/2, l=(0.5,0.5)) the cross-moment
    mass and bulk functionals are flat yet lambda_max(M_aug) separates
    UPWARD - falsifies the probe-in-disguise guess."""
    th = np.pi / 2
    rng = np.random.default_rng(23)
    q0, q1, mx0, mx1 = [], [], [], []
    for _ in range(22):
        for g, ql, ml in ((0.0, q0, mx0), (3.2, q1, mx1)):
            X, Y = _draw([0.5, 0.5], 0.2, 1600, th, rng, g)
            Xc = X - X.mean(0, keepdims=True)
            Yc = Y - Y.mean()
            Ys = Yc / np.std(Yc)
            from simulator import spectrum

            d, V = spectrum(Xc)
            b = Xc.T @ Ys / X.shape[0]
            ql.append(float(np.sum((V.T @ b) ** 2)))
            mx, _ = _aug_stats(X, Y)
            ml.append(mx)
    assert abs(np.mean(q0) - np.mean(q1)) < 0.06, (np.mean(q0), np.mean(q1))
    auc_up = _auc_greater(np.array(mx0), np.array(mx1))
    assert auc_up > 0.88, auc_up


def test_FC_low_edge_channel_c5_sub():
    """F-C / prediction P-1: at c=5 sub the BOTTOM edge separates at g=1
    (coupling pushes the low outlier down)."""
    l = [0.5 * np.sqrt(5.0)] * 3
    _, mns = _mx_cell(l, 5.0, 500, np.pi / 6, [0.0, 1.0], reps=45, seed=13)
    auc_bot_down = 1.0 - _auc_greater(mns[0.0], mns[1.0])
    assert auc_bot_down > 0.75, auc_bot_down


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
