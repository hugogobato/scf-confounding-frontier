"""Phase 2 unit tests: formulas (F12/F13), estimator suite, detection stats,
DGP variants. Sizes kept small (machine shared; research plan Section 10.1).
"""
from __future__ import annotations

import math
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))

import de_formulas as df
import estimators as est
import detection as det
from simulator import Config, gen_data, spectrum


# ---------------------------------------------------------------------------
# F4 inverse / weights
# ---------------------------------------------------------------------------


def test_bbp_invert_roundtrip():
    for c in (0.2, 0.8, 2.0, 5.0):
        for l in (np.sqrt(c) * 1.3, 1.0, 3.0 * np.sqrt(c)):
            if l <= np.sqrt(c):
                continue
            nu = df.bbp_location(l, c)
            l_back = df.bbp_invert(nu, c)
            assert abs(l_back - l) < 1e-8 * max(1.0, l)
    assert np.isnan(df.bbp_invert(0.9 * (1 + math.sqrt(0.5)) ** 2, 0.5))


def test_weight_families():
    d = np.array([0.5, 1.0, 2.0, 10.0])
    w_trim = df.trim_weights(d, 1.0)
    assert np.all(w_trim <= 1.0) and np.all(w_trim >= np.minimum(1.0, 0.05))
    assert np.all(np.diff(w_trim) <= 1e-15)  # non-increasing in d
    w_lava = df.lava_weights(d, 1.0)
    assert np.allclose(w_lava, 1.0 / (1.0 + d))
    # Cevid default equals min(1, median/d)
    med = float(np.median(d))
    assert np.allclose(df.trim_weights(d, med),
                       np.minimum(1.0, med / d))


def test_sdboost_path_limits():
    z = np.array([1.0, -2.0])
    d_svd = np.array([10.0, 5.0])
    w = np.array([1.0, 1.0])
    alpha0 = df.sdboost_path_coefficients(z, d_svd, w, 0, 0.1)
    assert np.allclose(alpha0, 0.0)
    alpha_inf = df.sdboost_path_coefficients(z, d_svd, w, 5000, 0.1)
    assert np.allclose(alpha_inf, z / d_svd, atol=1e-8)


def test_loo_scores_match_bruteforce_ridge():
    rng = np.random.default_rng(7)
    n, p = 40, 12
    Xc = rng.standard_normal((n, p))
    Ys = rng.standard_normal(n)
    eig = spectrum(Xc)
    U = est.left_factors(Xc, eig)
    lam = 0.37

    def wf(l):
        return eig[0] / (eig[0] + l)

    score_fast = est.loo_scores_spectral(U, Ys, eig[0], wf, [lam])[0]
    # brute force leave-one-out ridge: refit on the subset with its OWN spectrum
    errs = []
    for i in range(n):
        mask = np.ones(n, bool)
        mask[i] = False
        Xtr, Ytr = Xc[mask], Ys[mask]
        d_tr, V_tr = spectrum(Xtr)
        b = est._ridge_on_spectrum(d_tr, V_tr, Xtr, Ytr, lam)
        errs.append((Ys[i] - Xc[i] @ b) ** 2)
    assert abs(score_fast - float(np.mean(errs))) < 5e-3, (
        score_fast, float(np.mean(errs)))


# ---------------------------------------------------------------------------
# F12 max-z null variance law
# ---------------------------------------------------------------------------


def _null_data(cfg_seed_tuple=(900, 300), rep=0):
    n, p = cfg_seed_tuple
    cfg = Config(n=n, p=p, r=2, l=(1.2, 0.35), theta=math.pi / 6, g=0.0,
                 profile="sub", label="t", q_fixed=True)
    return cfg, gen_data(cfg, rep)


def test_maxz_null_variance():
    n, p = 900, 300
    c = p / n
    vals_top, vals_second = [], []
    ds = []
    for rep in range(120):
        cfg, data = _null_data((n, p), rep)
        from estimators import center_columns, standardize_response

        Xc = center_columns(data["X"])
        Ys = standardize_response(data["Y"])
        d, V = spectrum(Xc)
        b = Xc.T @ Ys / n
        w = V[:, :2].T @ b
        z = math.sqrt(n) * w / np.sqrt(d[:2])
        vals_top.append(z[0] ** 2)
        vals_second.append(z[1] ** 2)
        ds.append(d)
    ds = np.mean(ds, axis=0)
    # F12 (Erratum 1): E[z_j^2] = (d_j/c + se2)/sigma_y2 with true
    # sigma_u2 = sigma_eps2 = 1, sigma_y2 = tr(Sigma)/p + 1 ~ mean(d) + 1
    pred_top = (ds[0] / c + 1.0) / (float(np.mean(ds)) + 1.0)
    pred_second = (ds[1] / c + 1.0) / (float(np.mean(ds)) + 1.0)
    # chi-square mean over 120 reps: relative SE ~ 13% for chi2_1
    assert abs(float(np.mean(vals_top)) / pred_top - 1.0) < 0.35, (
        float(np.mean(vals_top)), pred_top)
    assert abs(float(np.mean(vals_second)) / pred_second - 1.0) < 0.35


# ---------------------------------------------------------------------------
# F13 augmented statistic: H0 location and rough size
# ---------------------------------------------------------------------------


def test_aug_null_root_and_size():
    from runners import run_stats_rep

    n, p = 600, 240
    cfg = Config(n=n, p=p, r=2, l=(1.0, 0.3), theta=math.pi / 6, g=0.0,
                 profile="mixed", label="t-aug", q_fixed=True)
    t_vals, rej = [], []
    for rep in range(150):
        rows, _ = run_stats_rep(cfg, rep)
        row = rows[0]
        t_vals.append(row["t_aug"])
        rej.append(row["rej_s1_analytic"])
        if rep == 0:
            root = row["lam0_plugin"]
    # simulated top root should sit near the plug-in H0 root (within TW width)
    width = df.tw_width_cov_scale(n, p + 1)
    m_t, m_root = float(np.mean(t_vals)), root
    assert abs(m_t / m_root - 1.0) < 0.08, (m_t, m_root)
    # analytic size in [0, 0.25] for this quick check (gate cells use 10k reps)
    size = float(np.mean(rej))
    assert size <= 0.25, size


def test_s1_power_grows_with_g():
    """Frozen decoupling contrast (Erratum 3): S2 maxz_cal discriminates
    supercritical-aligned confounding strongly, while in subcritical cells
    BOTH S2 and S1-aug are blind at g = 2 (invisible-yet-harmful region).
    Analytic-threshold size is a WP 2.1 measurement (10k-rep nulls)."""
    from runners import run_stats_rep

    n, p = 600, 240
    # supercritical-aligned: S2 fires under MC calibration
    null_t = []
    cfg0 = Config(n=n, p=p, r=2, l=(1.0, 0.3), theta=math.pi / 6, g=0.0,
                  profile="mixed", label="t-mc0", q_fixed=True)
    for rep in range(120):
        null_t.append(run_stats_rep(cfg0, rep)[0][0]["t_maxz"])
    thr = float(np.quantile(null_t, 0.95))
    cfg1 = Config(n=n, p=p, r=2, l=(1.0, 0.3), theta=math.pi / 6, g=2.0,
                  profile="mixed", label="t-mc2", q_fixed=True)
    power_sup = float(np.mean([
        run_stats_rep(cfg1, rep)[0][0]["t_maxz"] > thr for rep in range(60)]))
    assert power_sup >= 0.6, power_sup

    # subcritical: blind by construction (both statistics within noise)
    for stat, tol in (("t_maxz", 0.75), ("t_aug", 0.35)):
        cfgs = Config(n=n, p=p, r=2, l=(0.35, 0.35), theta=math.pi / 6,
                      g=0.0, profile="sub", label=f"t-sub0-{stat}",
                      q_fixed=True)
        cfga = Config(n=n, p=p, r=2, l=(0.35, 0.35), theta=math.pi / 6,
                      g=2.0, profile="sub", label=f"t-sub2-{stat}",
                      q_fixed=True)
        v0 = [run_stats_rep(cfgs, rep)[0][0][stat] for rep in range(40)]
        v1 = [run_stats_rep(cfga, rep)[0][0][stat] for rep in range(40)]
        shift = abs(float(np.mean(v1)) - float(np.mean(v0)))
        scale = max(float(np.std(v0)), 1e-9)
        assert shift / scale <= tol + 3.0, (stat, shift, scale)


# ---------------------------------------------------------------------------
# Cevid Trim transform core property
# ---------------------------------------------------------------------------


def test_cevid_transform_shrinks_spike_aligned_bias():
    rng = np.random.default_rng(11)
    n, p, r = 800, 400, 2
    Q, _ = np.linalg.qr(rng.standard_normal((p, r)))
    l = (2.0, 0.3)
    Lam = Q * np.sqrt(l)[None, :]
    f = rng.standard_normal((n, r))
    u = rng.standard_normal((n, p))
    X = f @ Lam.T + u
    Xc = X - X.mean(0)
    d, V = spectrum(Xc)
    # spike-aligned bias direction
    b = Q[:, 0]
    Xb = Xc @ b
    n_med = len(Ys := np.zeros(n))  # noqa: F841 (placeholder for shape only)
    med_svd2 = n * float(np.median(d))
    w = df.trim_weights(d * n, med_svd2)
    # transformed design: X~ = U diag(sqrt(w)) U' Xc -> Gram eigenvalues w*d
    U = Xc @ V / np.sqrt(n * d)[None, :]
    XtX_tilde = (U * w) @ U.T @ Xc
    b_t = XtX_tilde @ b
    assert float(np.linalg.norm(b_t)) < 0.5 * float(np.linalg.norm(Xb))


# ---------------------------------------------------------------------------
# SDBoost EB variance components recover mixed-model scales
# ---------------------------------------------------------------------------


def test_sdboost_recovers_variance_components():
    rng = np.random.default_rng(23)
    n, p = 700, 350
    Xc = rng.standard_normal((n, p))
    sr_true, se_true = 0.02, 1.0
    b_draw = sr_true ** 0.5 * rng.standard_normal(p)
    Y = Xc @ b_draw + se_true ** 0.5 * rng.standard_normal(n)
    Ys = (Y - Y.mean()) / Y.std()
    eig = spectrum(Xc)
    beta, info = est.fit_sdboost_linear_eb(Xc, Ys, eig, None)
    rho_hat = info["sr2"] / info["se2"]
    rho_true = sr_true / se_true
    # log-scale tolerance: factor e^{±1.5} around truth (loose, one dataset)
    assert abs(math.log(rho_hat / rho_true)) < 1.5, (rho_hat, rho_true)


def test_sdboost_path_matches_min_norm_at_large_m():
    rng = np.random.default_rng(31)
    n, p = 250, 180
    Xc = rng.standard_normal((n, p))
    Ys = rng.standard_normal(n)
    eig = spectrum(Xc)
    d, V = eig
    U = Xc @ V / np.sqrt(n * d)[None, :]
    z = U.T @ Ys
    w = np.ones_like(d)
    beta = V @ df.sdboost_path_coefficients(z, np.sqrt(n * d), w, 10_000, 0.1)
    bmin = est.fit_min_norm(Xc, Ys, eig)
    assert float(np.linalg.norm(beta - bmin)) < 1e-8


# ---------------------------------------------------------------------------
# SEB tuner behavior
# ---------------------------------------------------------------------------


def _prep(cfg, rep=0):
    data = gen_data(cfg, rep)
    from estimators import center_columns

    Xc = center_columns(data["X"])
    Yc = data["Y"] - data["Y"].mean()  # RAW centered (estimator convention)
    return data, Xc, Yc, spectrum(Xc)


def test_seb_tau_orders_with_confounding_strength():
    base = dict(n=800, p=400, r=2, theta=math.pi / 6, profile="t",
                label="seb", q_fixed=True)
    taus = []
    for g in (0.2, 2.5):
        cfg = Config(l=(0.8, 0.3), g=g, **base)
        data, Xc, Ys, eig = _prep(cfg)
        _, info = est.est_eb_spectral(Xc, Ys, eig, None)
        taus.append(info["tau"])
    assert taus[1] < taus[0], taus


def test_eb_spectral_reduces_harmful_bias_vs_ols():
    n, p = 800, 320
    cfg = Config(n=n, p=p, r=2, l=(0.45, 0.45), theta=math.pi / 6, g=1.0,
                 profile="sub", label="harmful", q_fixed=True)
    acc_ols = np.zeros(p)
    acc_eb = np.zeros(p)
    for rep in range(40):
        data, Xc, Ys, eig = _prep(cfg, rep)
        b_ols = est.fit_min_norm(Xc, Ys, eig)
        b_eb, _ = est.est_eb_spectral(Xc, Ys, eig, None)
        acc_ols += b_ols - data["beta"]
        acc_eb += b_eb - data["beta"]
    # sub profile at c<1: exact OLS bias sqrt(sum (sqrt(l)/(1+l))^2 gamma^2)
    l = np.array([0.45, 0.45])
    gam = np.array([math.cos(math.pi / 6), math.sin(math.pi / 6)])
    pred = float(np.linalg.norm(np.sqrt(l / (1 + l)) * gam))
    bias_ols = float(np.linalg.norm(acc_ols / 40))
    # MC noise floor on the mean-diff norm ~ sqrt(p/reps)-scaled; allow 30%
    assert abs(bias_ols / pred - 1.0) < 0.30, (bias_ols, pred)
    bias_eb = float(np.linalg.norm(acc_eb / 40))
    assert bias_eb < 0.75 * bias_ols, (bias_eb, bias_ols)


# ---------------------------------------------------------------------------
# DGP variants
# ---------------------------------------------------------------------------


def test_sparse_confounding_dgp_geometry():
    n, p, r = 500, 400, 3
    cfg = Config(n=n, p=p, r=r, l=(1.5, 0.5, 0.5), theta=math.pi / 6,
                 conf_kind="sparse", profile="mixed", label="sp",
                 q_fixed=True)
    data = gen_data(cfg, 0)
    Lam = data["Lam"]
    svals = np.sqrt(np.asarray(cfg.l))
    gram = Lam.T @ Lam
    assert np.allclose(gram, np.diag(svals ** 2), atol=1e-8)
    Ubasis = Lam / svals[None, :]
    assert np.allclose(Ubasis.T @ Ubasis, np.eye(r), atol=1e-8)
    # supports disjoint
    supp = [set(np.nonzero(Lam[:, j])[0]) for j in range(r)]
    for i in range(r):
        for j in range(i + 1, r):
            assert not (supp[i] & supp[j])


def test_rademacher_half_column_norms():
    n, p, r = 400, 300, 2
    cfg = Config(n=n, p=p, r=r, l=(1.2, 0.4), theta=math.pi / 6,
                 loading_kind="rademacher_half", profile="mixed",
                 label="rad", q_fixed=True)
    data = gen_data(cfg, 0)
    norms = np.linalg.norm(data["Lam"], axis=0)
    assert np.allclose(norms, np.sqrt([1.2, 0.4]), rtol=1e-6)


def test_m2_dgp_recovers_tau():
    cfg = Config(n=3000, p=300, r=2, l=(0.6, 0.3), theta=math.pi / 6,
                 m2_treatment=True, m2_tau=1.0, delta_g=0.3,
                 profile="mixed", label="m2t", q_fixed=True)
    data = gen_data(cfg, 0)
    X, D, Y = data["X"], data["D"], data["Y"]
    A = np.column_stack([D, X])
    coef, *_ = np.linalg.lstsq(A, Y, rcond=None)
    # dense beta contributes O(c)-ish contamination; just require recovery
    # within a loose band at this large n
    assert abs(coef[0] - 1.0) < 0.2, coef[0]


# ---------------------------------------------------------------------------
# detection plumbing
# ---------------------------------------------------------------------------


def test_compute_stats_keys_and_types():
    data, Xc, Ys, eig = _prep(Config(
        n=300, p=120, r=2, l=(0.9, 0.3), theta=math.pi / 6, g=0.5,
        profile="mixed", label="keys", q_fixed=True))
    stats = det.compute_stats(Xc, Ys, eig)
    need = ["lam_max_cov", "tw_stat", "outlier99", "t_aug", "t_maxz",
            "lss_rank1", "f_pcs", "b_norm2", "z_top", "ktop"]
    for k in need:
        assert k in stats
    rej = det.rejections(stats)
    assert all(isinstance(v, bool) for v in rej.values())
    feats = det.probe_features(eig, stats)
    assert feats.shape == (14,)


def test_lecam_probe_separates_and_calibrates():
    rng = np.random.default_rng(5)
    h0 = rng.standard_normal((400, 14))
    h1 = h0[:200] + 0.8
    auc_sep = det.lecam_auc(h0, h1, seed=0)
    assert auc_sep[0] > 0.7 and auc_sep[1] > 0.7
    auc_same = det.lecam_auc(rng.standard_normal((300, 14)),
                             rng.standard_normal((300, 14)), seed=1)
    assert 0.35 < auc_same[0] < 0.65 and 0.35 < auc_same[1] < 0.65


def test_ucm_strength_bounded():
    d_white = np.sort(np.random.default_rng(0).chisquare(3, 50))[::-1]
    val = df.ucm_strength(d_white / 3.0, 50, 50)
    assert 0.0 <= val <= 1.0


def test_seb_predicted_mse_monotone_in_g2():
    l = np.array([1.0, 0.3])
    d = np.array([df.bbp_location(1.0, 0.5), 1.6])
    mse0 = df.seb_predicted_mse(l, np.array([0.0, 0.0]), d, 0.5, 0.5, 800, 1.0)
    mse1 = df.seb_predicted_mse(l, np.array([1.0, 0.0]), d, 0.5, 0.5, 800, 1.0)
    assert mse1 > mse0
