"""WP 1.4 identity tests: every transcribed formula gets a numerical self-check.

Deviations from the plan's literal tolerances are documented in
docs/de_formula_sheet.md and docs/pilot_memo.md: RMT simulation checks run at
n <= 8000 (memory guard, research plan Section 10.1) instead of n = 20000,
which is compensated by multi-rep averaging and a convergence-trend assertion.
"""
from __future__ import annotations

import math
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.sparse.linalg import LinearOperator, eigsh

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))

import de_formulas as df
from simulator import Config, fit_ols, fit_ridge, gen_data, run_rep, spectrum


# --------------------------------------------------------------------------
# F1/F2/F3: exact population identities (algebraic)
# --------------------------------------------------------------------------


def _random_setup(seed=0, p=60, r=4):
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.standard_normal((p, r)))
    l = np.sort(rng.uniform(0.1, 3.0, size=r))[::-1]
    gamma = rng.standard_normal(r)
    beta = rng.standard_normal(p)
    return Q, l, gamma, beta, rng


def test_population_ols_identity():
    Q, l, gamma, beta, rng = _random_setup()
    sigma2 = 1.7
    Lam = Q * (math.sqrt(sigma2) * np.sqrt(l))
    Sigma = sigma2 * np.eye(len(beta)) + Lam @ Lam.T
    cov_xy = Sigma @ beta + Lam @ gamma
    # identity 1 of plan Section 2.3: plim beta_OLS - beta = Sigma^{-1} Lambda gamma
    lhs = np.linalg.solve(Sigma, cov_xy) - beta
    rhs = np.linalg.solve(Sigma, Lam @ gamma)
    assert np.allclose(lhs, rhs, atol=1e-10)
    # component form matches the matrix form, and de_formulas agrees
    comp = sum(
        math.sqrt(l[j] / sigma2) / (1.0 + l[j]) * gamma[j] * Q[:, j]
        for j in range(len(l))
    )
    assert np.allclose(comp, rhs, atol=1e-10)
    # de_formulas returns the r coefficients on (u_1..u_r); mapping through Q
    # must reproduce the p-vector
    coeffs = df.ols_bias_vector(l, gamma, sigma2)
    assert np.linalg.norm(Q @ coeffs - comp) < 1e-10


def test_population_ridge_identity():
    Q, l, gamma, beta, rng = _random_setup(seed=1)
    sigma2 = 1.3
    Lam = Q * (math.sqrt(sigma2) * np.sqrt(l))
    Sigma = sigma2 * np.eye(Q.shape[0]) + Lam @ Lam.T
    for lam in (0.05, 0.5, 5.0):
        pred = np.linalg.solve(Sigma + lam * np.eye(Q.shape[0]), Lam @ gamma)
        got = df.ridge_bias_vector(l, gamma, lam, sigma2)
        # component form: sum_j sqrt(l_j sigma2)/(sigma2(1+l_j)+lam) gamma_j u_j
        comp = sum(math.sqrt(l[j] * sigma2) / (sigma2 * (1 + l[j]) + lam) * gamma[j]
                   * Q[:, j] for j in range(len(l)))
        assert np.allclose(pred, comp, atol=1e-10)
        assert np.linalg.norm(Q @ got - comp) < 1e-10


def test_pca_trim_tradeoff():
    Q, l, gamma, beta, _ = _random_setup(seed=2, p=50, r=4)
    sigma2 = 1.0
    k = 2
    bias = df.pca_trim_bias_vector(l, gamma, k, sigma2)
    expected = sum(math.sqrt(l[j] * sigma2) / (sigma2 * (1 + l[j])) * gamma[j]
                   * Q[:, j] for j in range(k))
    assert np.linalg.norm(Q @ bias - expected) < 1e-10
    # trimming removes exactly the dropped-direction bias
    dropped = sum(sigma2 * l[j] / (sigma2 * (1 + l[j])) * gamma[j] * Q[:, j]
                  for j in range(k, len(l)))
    assert np.linalg.norm(dropped) > 1e-6


def test_minnorm_isotropic_limit():
    # F8 anchor: r = 0 reduces to E[beta_hat] = beta/c
    for c in (2.0, 5.0):
        got = df.minnorm_total_bias_norm(np.array([]), np.array([]), c, p=1000)
        assert abs(got - abs(1.0 / c - 1.0)) < 1e-10


def test_ridge_capture_reduces_to_superseded_capture():
    # UPDATED 2026-08-25: the proved shifted-resolvent ridge capture
    # (docs/theory_T1_capture_law.md Section 6) has lam -> 0 limit equal to
    # the PROVED min-norm capture law (1+l)/(c+l), NOT the superseded
    # xi-based guess. The old assertion here enshrined the xi-split
    # interpolation, which the decisive reconciliation check falsified
    # (max err 0.081 vs 0.003 for the proved form); the superseded form is
    # preserved as ridge_capture_superseded and still reduces to
    # bgn_capture_superseded as lam -> 0.
    l = np.array([2.0, 1.5])
    c = 5.0
    cap_lam = df.ridge_capture(l, 1e-12, c)
    assert np.allclose(cap_lam, df.minnorm_capture(l, c), rtol=1e-6)
    cap_old = df.ridge_capture_superseded(l, 1e-12, c)
    assert np.allclose(cap_old, df.bgn_capture_superseded(l, c), rtol=1e-6)


def test_minnorm_capture_law_anchors():
    # pilot-validated conjecture cap = (1 + l)/(c + l): boundary behavior
    l = np.array([0.0, 3.0])
    c = 5.0
    cap = df.minnorm_capture(l, c)
    assert abs(cap[0] - 1 / c) < 1e-12          # l -> 0: uniform rowspace
    assert np.allclose(cap, (1 + l) / (c + l), rtol=1e-14)
    assert df.minnorm_capture(np.array([1e9]), c)[0] > 1 - 1e-6  # l -> inf


def test_mp_stieltjes_inv_matches_quadratic():
    for lam in (0.01, 0.1, 1.0, 10.0):
        m = df.mp_stieltjes_inv(lam, 0.5)
        resid = lam * m**2 + (lam + 0.5) * m - 1.0
        assert abs(resid) < 1e-10


# --------------------------------------------------------------------------
# F4/F5/F7: classical laws checked against simulation
# --------------------------------------------------------------------------


def _top_eigs_spiked(n: int, p: int, l: float, reps: int, seed=0):
    """Sample lambda_max (and top eigenvector overlap) for one spike l.

    Uses orthogonally-invariant reduction X = sqrt(l) z q' + W with q = e_1
    (i.e., only column 0 carries the spike), which has exactly the same
    eigenvalue and eigenvector-overlap laws as any Haar loading direction.
    """
    rng = np.random.default_rng(seed)
    lams, ovls = [], []
    for _ in range(reps):
        z = rng.standard_normal(n)
        X = rng.standard_normal((n, p))
        X[:, 0] += math.sqrt(l) * z
        Aop = LinearOperator(
            (p, p), matvec=lambda v: X.T @ (X @ v) / n, dtype=float
        )
        vals, vecs = eigsh(Aop, k=1, which="LA", tol=1e-9, maxiter=5000,
                           v0=np.ones(p) / math.sqrt(p))
        lams.append(vals[0])  # Aop already carries the 1/n normalization
        q = np.zeros(p)
        q[0] = 1.0
        ovls.append(float((vecs[:, 0] @ q) ** 2))
    return np.mean(lams), np.mean(ovls), np.std(lams)


def test_bbp_location():
    c, l = 0.5, 2.0
    mu = df.bbp_location(l, c)
    assert abs(mu - 3.75) < 1e-12  # corrected formula, not the plan's 3.0 slip
    m_small, _, s_small = _top_eigs_spiked(4000, 2000, l, reps=12, seed=1)
    m_big, _, s_big = _top_eigs_spiked(8000, 4000, l, reps=12, seed=2)
    rel_small = abs(m_small - mu) / mu
    rel_big = abs(m_big - mu) / mu
    # finite-size tolerance 4 percent + convergence trend toward mu
    assert rel_small < 0.04, f"BBP location off at n=4000: {rel_small:.3%}"
    assert rel_big < 0.04, f"BBP location off at n=8000: {rel_big:.3%}"
    mc_se_big = s_big / math.sqrt(12) / mu
    assert rel_big < max(0.04, 4 * mc_se_big)


def test_bgn_overlap():
    c, l = 0.5, 2.0
    xi = df.bgn_overlap(l, c)
    assert abs(xi - 0.7) < 1e-12
    _, ovl_mean, _ = _top_eigs_spiked(6000, 3000, l, reps=12, seed=3)
    assert abs(ovl_mean - xi) < 0.02, f"overlap {ovl_mean:.4f} vs xi {xi}"


def test_subcritical_no_outlier():
    # subcritical spike sticks near the bulk edge, overlap ~ 0
    c, l = 0.5, 0.4
    assert df.bbp_location(l, c) == (1 + math.sqrt(c)) ** 2
    assert df.bgn_overlap(l, c) == 0.0
    m, ovl, _ = _top_eigs_spiked(3000, 1500, l, reps=8, seed=4)
    edge = (1 + math.sqrt(c)) ** 2
    assert abs(m - edge) / edge < 0.03
    assert ovl < 0.06


def test_tw_calibration():
    n, p = 400, 200
    thr95 = df.tw_threshold(n, p, 1.0, q=df.TW1_Q95)
    thr99 = df.tw_threshold(n, p, 1.0, q=df.TW1_Q99)
    rng = np.random.default_rng(11)
    rej95 = rej99 = 0
    reps = 5000
    for _ in range(reps):
        X = rng.standard_normal((n, p))
        lam_max = np.linalg.eigvalsh(X.T @ X / n)[-1]
        rej95 += lam_max > thr95
        rej99 += lam_max > thr99
    rate95, rate99 = rej95 / reps, rej99 / reps
    assert 0.04 <= rate95 <= 0.06, f"TW q95 calibration off: {rate95}"
    assert 0.005 <= rate99 <= 0.02, f"TW q99 calibration off: {rate99}"


# --------------------------------------------------------------------------
# Simulator wiring checks
# --------------------------------------------------------------------------


def test_ols_exact_finite_n():
    """F1 exact identity (conditional on Q): mean(beta_hat - beta) = a."""
    cfg = Config(n=1200, p=300, r=2, l=(1.5, 0.5), theta=math.pi / 6,
                 q_fixed=True, profile="check", label="exact")
    gam = np.array([math.cos(cfg.theta), math.sin(cfg.theta)])
    pred = df.ols_bias_vector(np.asarray(cfg.l), gam, 1.0)
    acc = np.zeros(cfg.p)
    reps = 1000
    Q_ref = None
    for rep in range(reps):
        data = gen_data(cfg, rep)
        Q_ref = data["Q"]
        d, V = spectrum(data["X"])
        bh = fit_ols(data["X"], data["Y"], (d, V))
        acc += bh - data["beta"]
    got = acc / reps
    pred_full = Q_ref @ pred  # factor coordinates -> p-vector
    pred_dir = pred_full / np.linalg.norm(pred_full)
    # sharp check: component along the predicted direction
    proj_rel = abs(float(got @ pred_dir) - np.linalg.norm(pred_full)) / \
        np.linalg.norm(pred_full)
    assert proj_rel < 0.05, f"projected exact-bias mismatch: {proj_rel:.3%}"
    # coarser check: total vector within MC noise of the prediction
    rel = np.linalg.norm(got - pred_full) / np.linalg.norm(pred_full)
    assert rel < 0.12, f"exact finite-n OLS bias mismatch: {rel:.3%}"


def test_ridge_zero_lambda_matches_ols():
    cfg = Config(n=300, p=120, r=2, l=(1.0, 0.4), theta=0.3,
                 profile="check", label="wiring")
    data = gen_data(cfg, 7)
    eig = spectrum(data["X"])
    assert np.linalg.norm(fit_ridge(data["X"], data["Y"], eig, 1e-12)
                          - fit_ols(data["X"], data["Y"], eig)) < 1e-7


def test_run_rep_schema():
    cfg = Config(n=250, p=100, r=2, l=(1.2, 0.5), theta=0.4, twin_gamma0=True,
                 profile="smoke", label="schema")
    rows, diffs = run_rep(cfg, 3, (0.1, 1.0))
    df_rows = pd.DataFrame(rows)
    for col in ("config_id", "rep", "estimator", "lambda", "k_select",
                "rel_err", "runtime_s", "lam_max_cov", "tw_stat", "outlier99",
                "bbp_pred", "xi1_pred", "overlap1", "overlap2"):
        assert col in df_rows.columns, col
    ests = set(df_rows["estimator"])
    assert {"ols", "ridge_fixed", "pca_onatski", "ols_g0",
            "ridge_fixed_g0"} <= ests
    assert df_rows["rel_err"].notna().all()
    assert np.isfinite(df_rows["tw_stat"]).all()
    # twins share everything but the gamma link, so their difference isolates it
    key = "ols"
    assert key in diffs and key + "_g0" in diffs
    assert np.isfinite(diffs[key]).all()


def test_ledger_hash_deterministic():
    docs = Path(__file__).resolve().parents[1] / "docs"
    h1 = df.ledger_hash(docs)
    h2 = df.ledger_hash(docs)
    assert h1 == h2 and len(h1) == 12
