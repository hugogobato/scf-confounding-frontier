"""T7 trimmed-tau tests (closed 2026-08-25, docs/theory_T7_trimmed_tau.md).

Falsifiers:
  - SM identity tau_hat_ols = d'G^-1 Y / (1 + d'G^-1 D)  vs pinv joint OLS
  - OLS shrinkage plim -tau/(1+Lambda_D) with the three-line Lambda_D DE
  - trim flatness across c (systematic channel scale, loose tol)
"""

import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "code"))


def m2_draw(l, c, n, rng):
    """One M2 draw matching model card / runner conventions."""
    p = int(round(c * n))
    r = len(l)
    Q, _ = np.linalg.qr(rng.standard_normal((p, r)))
    gam = np.zeros(r)
    gam[0] = np.cos(np.pi / 6)
    if r > 1:
        gam[1] = np.sin(np.pi / 6)
    sv = np.sqrt(np.asarray(l))
    beta = rng.standard_normal(p)
    beta /= np.linalg.norm(beta)
    f = rng.standard_normal((n, r))
    U = rng.standard_normal((n, p))
    eps = rng.standard_normal(n)
    X = (f * sv[None, :]) @ Q.T + U
    kpi = max(3, p // 100)
    pi = np.zeros(p)
    pi[:kpi] = 1.0 / np.sqrt(kpi)
    delta = rng.standard_normal(r)
    delta *= 0.3 / max(float(np.linalg.norm(delta)), 1e-12)
    nu = rng.standard_normal(n)
    D = X @ pi + f @ delta + nu
    Y = D + X @ beta + f @ gam + eps
    return X, D, Y, pi, delta


def test_ols_sm_identity():
    """tau_hat min-norm joint OLS equals d'G^-1Y/(1+d'G^-1D) exactly."""
    rng = np.random.default_rng(7)
    l = [4.243, 0.707, 0.707, 0.707, 0.707]
    n, c = 250, 2.0
    X, D, Y, _, _ = m2_draw(l, c, n, rng)
    Xc = X - X.mean(0, keepdims=True)
    dc = D - D.mean()
    yc = Y - Y.mean()
    M = np.column_stack([dc, Xc])
    a_pinv = float((np.linalg.pinv(M) @ yc)[0])
    G = Xc @ Xc.T
    num = float(dc @ np.linalg.solve(G + np.outer(dc, dc), yc))
    assert abs(a_pinv - num) < 1e-8 * max(1.0, abs(num))


def test_ols_shrinkage_plim():
    """OLS tau bias ~ -tau/(1+Lambda_D) at c=2; Lambda_D three-line formula."""
    l = np.array([4.243, 0.707, 0.707, 0.707, 0.707])
    c, n, reps = 2.0, 400, 60
    t = 1.0 / (c - 1.0)
    lam_d_pred = (1.0 - 1.0 / c) + float(
        np.sum((t * (1.0 + t) / (1.0 + t * (1.0 + l))) * (0.3 ** 2 / len(l)))
    ) + 1.0
    rng = np.random.default_rng(11)
    errs, lam_ds = [], []
    for _ in range(reps):
        X, D, Y, _, _ = m2_draw(list(l), c, n, rng)
        Xc = X - X.mean(0, keepdims=True)
        dc = D - D.mean()
        yc = Y - Y.mean()
        G = Xc @ Xc.T
        one = np.ones(n) / n
        big = 1e6 * float(np.trace(G) / n)
        L = np.linalg.cholesky(G + big * np.outer(one, one))

        def Ginv(v):
            y = np.linalg.solve(L, v)
            return np.linalg.solve(L.T, y)

        num = float(dc @ Ginv(yc))
        den = 1.0 + float(dc @ Ginv(dc))
        errs.append(num / den - 1.0)
        lam_ds.append(den - 1.0)
    lam_meas = float(np.mean(lam_ds))
    err_mean = float(np.mean(errs))
    assert abs(lam_meas - lam_d_pred) < 0.25, (lam_meas, lam_d_pred)
    assert abs(err_mean - (-1.0 / (1.0 + lam_d_pred))) < 0.12, (
        err_mean,
        -1.0 / (1.0 + lam_d_pred),
    )


def test_trim_flat_across_c():
    """Trim systematic-error SD is flat in c while OLS error grows."""
    rng = np.random.default_rng(13)
    out = {}
    for c in (0.5, 2.0):
        l = [3.0 * np.sqrt(c)] + [0.5 * np.sqrt(c)] * 4
        n = 300 if c < 1 else 450
        p = int(round(c * n))
        errs_trim, errs_ols = [], []
        for _ in range(40):
            X, D, Y, _, _ = m2_draw(l, c, n, rng)
            Xc = X - X.mean(0, keepdims=True)
            dc = D - D.mean()
            yc = Y - Y.mean()
            d_spec = np.linalg.eigvalsh(
                Xc @ Xc.T / n if c < 1 else Xc @ Xc.T / n
            )
            from simulator import spectrum

            d_srt, V = spectrum(Xc)
            V = V[:, :5]
            S = Xc @ V
            Ms = np.column_stack([dc, S])
            errs_trim.append(float((np.linalg.lstsq(Ms, yc, rcond=None)[0])[0]) - 1.0)
            Mo = np.column_stack([dc, Xc])
            errs_ols.append(float((np.linalg.pinv(Mo) @ yc)[0]) - 1.0)
        out[c] = (np.std(errs_trim), np.abs(errs_ols).mean() if False
                  else float(np.mean(np.abs(errs_ols))))
    sd_sub, sd_sup = out[0.5][0], out[2.0][0]
    assert abs(sd_sup - sd_sub) < 0.6 * sd_sub, (sd_sub, sd_sup)  # flatness
    assert out[2.0][1] > 2.5 * sd_sup  # OLS inflation at c=2


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
