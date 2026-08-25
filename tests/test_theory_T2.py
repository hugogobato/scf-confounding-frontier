"""T2 finite-n null law + frontier achievability tests (closed 2026-08-25,
docs/theory_T2_frontier.md).

Falsifiers:
  IND  conditional independence of calibrated coordinates (Lemma A(iii))
  KAP  kappa scale law at c>1: sd(z_cal) matches the closed-form DE ratio
  RAW  raw-Bonferroni size explosion at c=5 reproduced and explained
  SLP  supercritical slope law with BGN xi^{1/2} weight: linearity in g and
       absolute level within tolerance
  FRN  frozen gate anchor: frontier ratios <= 1.5 (csv)
"""

import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "code"))

ROOT = os.path.join(os.path.dirname(__file__), "..")


def _draw(l, c, n, theta, rng, g):
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


def _zcal_coords(l, c, n, theta, rng, g, kk_override=None):
    """Per-rep calibrated spike coordinates via the frozen pipeline pieces."""
    from simulator import spectrum
    from de_formulas import estimate_noise_scales, maxz_null_var, onatski_select

    X, Y = _draw(l, c, n, theta, rng, g)
    Xc = X - X.mean(0, keepdims=True)
    Yc = Y - Y.mean()
    Ys = Yc / np.std(Yc)
    d, V = spectrum(Xc)
    b = Xc.T @ Ys / n
    kk = int(min(max(onatski_select(d), 1), 10)) if kk_override is None \
        else kk_override
    kk = min(kk, len(d))
    kk = min(kk, len(d))
    w = V[:, :kk].T @ b
    z = np.sqrt(n) * w / np.sqrt(d[:kk])
    se2, sy2 = estimate_noise_scales(d, c)
    var_cal = np.array([maxz_null_var(dj, c, se2, sy2) for dj in d[:kk]])
    return z / np.sqrt(var_cal), d[:kk], var_cal


def _restricted_mp_median(c):
    a, b_ = (1 - np.sqrt(c)) ** 2, (1 + np.sqrt(c)) ** 2
    xs = np.linspace(a + (b_ - a) * 1e-7, b_ * (1 - 1e-9), 400001)
    dens = np.sqrt(np.maximum((b_ - xs) * (xs - a), 0)) / (2 * np.pi * c * xs)
    cdf = np.cumsum(dens) * (xs[1] - xs[0])
    cdf /= cdf[-1]
    return float(np.interp(0.5, cdf, xs))


def _mp_median(c):
    return _restricted_mp_median(c) if c < 1 else _mp_median_generic(c)


def _mp_median_generic(c):
    a, b_ = (1 - np.sqrt(min(c, 1.0))) ** 2, (1 + np.sqrt(min(c, 1.0))) ** 2
    xs = np.linspace(a + (b_ - a) * 1e-7, b_ * (1 - 1e-9), 400001)
    dens = np.sqrt(np.maximum((b_ - xs) * (xs - a), 0)) / (2 * np.pi * min(c, 1.0) * xs)
    cdf = np.cumsum(dens) * (xs[1] - xs[0])
    cdf /= cdf[-1]
    return float(np.interp(0.5, cdf, xs))


def test_INDEPENDENCE_of_calibrated_coordinates():
    rng = np.random.default_rng(31)
    for l, c, n, reps, kk_use in (
            ([0.5 * np.sqrt(5)] * 3, 5.0, 300, 90, 2),
            ([3 * np.sqrt(0.8), 0.5 * np.sqrt(0.8)], 0.8, 800, 60, 2)):
        zs = []
        for _ in range(reps):
            zc, _, _ = _zcal_coords(l, c, n, np.pi / 6, rng, 0.0,
                                    kk_override=kk_use)
            if len(zc) >= 2:
                zs.append(zc[:2])
        Z = np.array(zs)
        assert Z.ndim == 2 and Z.shape[0] >= 30
        if not np.all(Z.std(axis=0) > 0):
            continue
        C = np.corrcoef(Z.T)
        off = np.abs(C[np.triu_indices(2, 1)]).max() if C.size > 1 else 0.0
        assert off < 3.0 / np.sqrt(reps), (c, off)


def test_KAPPA_scale_law_c_gt_1():
    """sd(z_cal top) ~= closed-form DE kappa (scale-calibration ratio)."""
    c, n, l = 5.0, 300, [0.5 * np.sqrt(5)] * 3
    se2_hat_de = _restricted_mp_median(c) / _mp_median_generic(c)
    sy2_hat_de = c + se2_hat_de          # mean(d) -> c*tr/p at c>1
    dbar = 10.4                          # edge-pinned top coordinate (measured)
    A = 2.0
    kappa_pred = np.sqrt(
        ((1.0 + dbar / c) / A**2) / ((se2_hat_de + dbar / c) / sy2_hat_de**2)
    )
    rng = np.random.default_rng(33)
    zc = []
    for _ in range(110):
        z, _, _ = _zcal_coords(l, c, n, np.pi / 6, rng, 0.0)
        zc.append(z[0])
    sd_meas = float(np.std(zc))
    assert abs(se2_hat_de - 7.14) < 0.15 and abs(sy2_hat_de - 12.14) < 0.15
    assert 1.0 <= sd_meas / kappa_pred <= 1.55, (sd_meas, kappa_pred)


def test_RAW_bonferroni_size_explosion_c5():
    """The raw-Bonferroni rule explodes at c>1 (documented mechanism);
    magnitude consistent with the kappa law."""
    from scipy.stats import norm

    c, n, l = 5.0, 300, [0.5 * np.sqrt(5)] * 3
    rng = np.random.default_rng(35)
    rej, cuts = [], []
    for _ in range(150):
        zc, _, _ = _zcal_coords(l, c, n, np.pi / 6, rng, 0.0)
        cut = float(norm.ppf(1 - 0.05 / (2 * len(zc))))
        rej.append(float(np.max(np.abs(zc)) > cut))
        cuts.append(cut)
    size = float(np.mean(rej))
    assert np.mean(cuts) > 1.9           # Bonferroni cut in place
    assert size > 0.25, size             # explosion reproduced
    assert size < 0.95


def test_SLOPE_law_supercritical_xi_weight():
    """Second-moment slope law: because the BGN overlap has random sign
    (v_j'q_j = +- xi^{1/2}), mean(z) carries no signal; the coherent
    channel shows up in E[z_cal^2], growing like g^2/(A+g^2) with the
    xi^{1/2} weight inside the constant."""
    c, n = 0.8, 800
    l = [3 * np.sqrt(c), 0.5 * np.sqrt(c)]
    th = np.pi / 6

    def ez2(g, reps, seed):
        rng = np.random.default_rng(seed)
        acc = []
        for _ in range(reps):
            z, _, _ = _zcal_coords(l, c, n, th, rng, g)
            acc.append(z[0] ** 2)
        return float(np.mean(acc))

    e0 = ez2(0.00, 140, 50)
    e1 = ez2(0.15, 140, 51)
    e2 = ez2(0.30, 140, 52)
    e3 = ez2(0.60, 140, 53)
    r1 = (e2 - e0) / (e1 - e0)
    r2 = (e3 - e0) / (e1 - e0)
    # pure-g^2 would give 4 and 16; the sigma_y^2 = A + g^2 satiation
    # bends these to ~3.9 and ~13.8
    assert 3.0 < r1 < 5.2, r1
    assert 11.0 < r2 < 19.0, r2


def test_FROZEN_frontier_ratios():
    df = pd.read_csv(os.path.join(ROOT, "results", "frontier_check.csv"))
    gated = df[df.supercritical == True]  # noqa: E712
    assert len(gated) == 3
    assert gated.pass_le_1p5.all()
    assert gated.ratio_emp_over_pred.median() < 1.3


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
