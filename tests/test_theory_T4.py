"""T4 hard-trim dominance tests (closed 2026-08-25,
docs/theory_T4_hard_trim_dominance.md).

Estimator class: spectral REGRESSION maps
    beta_hat_m = V diag(m_j(d_hat)) V' b_raw,   m_j >= 0,
i.e. members are min-norm OLS (m = 1/d), ridge (m = 1/(d+lam)), hard trims
(m = 1{j in A}/d), soft-trim variants (m = 1{j in A}/(d+s)), lava-shaped
(m = 1/(d(1+lam d))). NOTE: V diag(w) V' b WITHOUT the 1/d factor is a
projection-type map, NOT an estimator of beta; the 1/d bookkeeping is
load-bearing (this exact pitfall was caught by the first falsifier run).

Falsifiers:
  SUB  all-subcritical cell: trim(k=0) kills the confounding channel;
       OLS transmits sqrt(l)/(c+l)*gamma (cap law); ridge transmits
       cap(lam)*sqrt(l)/(1+l)*gamma (T1.c law); lava-shaped soft transmits
       strictly more than trim => dominance on the invisible cell
  SUP  mixed cell: retained supercritical direction passes the ideal F3
       value sqrt(l)/(1+l)*gamma; dropped subcritical direction ~0; OLS
       passes the cap-law value on both
  ANCH frozen-csv anchors of the OLS transmission (two cells, <3%)
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
from de_formulas import ridge_capture  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")


def _draw_cell(Q, sv, r, n, theta, rng, g):
    p = Q.shape[0]
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
    return X - X.mean(0, keepdims=True), Y - Y.mean()


def _spectral_maps(Xc, Yc):
    """Regression maps m for the roster; RAW centered response (D3)."""
    from simulator import spectrum

    d, V = spectrum(Xc)
    b = V.T @ (Xc.T @ Yc / Xc.shape[0])
    out = {
        "ols": V @ (b / d),
        "ridge05": V @ (b / (d + 0.5)),
        "lava": V @ (b / (d * (1.0 + 0.5 * d))),
        "trim0": np.zeros(Xc.shape[1]),
    }
    return out


def _twin_means(l, c, n, theta, reps, seed):
    rng = np.random.default_rng(seed)
    p = int(round(c * n))
    r = len(l)
    Q, _ = np.linalg.qr(rng.standard_normal((p, r)))
    sv = np.sqrt(np.asarray(l))
    keys = ("ols", "ridge05", "lava", "trim0")
    acc = {g: {k: np.zeros(r) for k in keys} for g in (0.0, 1.0)}
    for g in (0.0, 1.0):
        for _ in range(reps):
            Xc, Yc = _draw_cell(Q, sv, r, n, theta, rng, g)
            est = _spectral_maps(Xc, Yc)
            for k, bh in est.items():
                acc[g][k] += Q.T @ bh
    return {g: {k: v / reps for k, v in acc[g].items()} for g in acc}


def _delta(means):
    return {k: means[1.0][k] - means[0.0][k] for k in means[0.0]}


def test_SUB_subcell_transmission_and_dominance():
    l = [0.5 * np.sqrt(2.0)] * 2
    c, n, th, g = 2.0, 500, np.pi / 6, 1.0
    means = _twin_means(l, c, n, th, reps=180, seed=77)
    dl = _delta(means)
    dirv = np.array([np.cos(th), np.sin(th)])
    lv = np.asarray(l)
    a_gamma = np.sqrt(lv) / (1.0 + lv) * g * dirv
    # trim(k=0) attains the zero floor
    assert np.linalg.norm(dl["trim0"]) < 0.02, dl["trim0"]
    # OLS transmits the cap law sqrt(l)/(c+l)
    pred_ols = np.sqrt(lv) / (c + lv) * g * dirv
    assert np.linalg.norm(dl["ols"] - pred_ols) < 0.04, (dl["ols"], pred_ols)
    # ridge transmits cap(lam) via the T1.c law
    cap_lam = float(ridge_capture(np.array([lv[0]]), 0.5, c)[0])
    pred_rid = cap_lam * a_gamma
    assert np.linalg.norm(dl["ridge05"] - pred_rid) < 0.04, (
        dl["ridge05"], pred_rid)
    # dominance: every soft family strictly loses to the trim here
    for k in ("ols", "ridge05", "lava"):
        assert np.linalg.norm(dl[k]) > 4.0 * max(
            np.linalg.norm(dl["trim0"]), 0.004), (k, dl[k])


def test_SUP_mixed_oracle_profile():
    l = [3 * np.sqrt(2.0), 0.5 * np.sqrt(2.0)]
    c, n, th, g = 2.0, 700, np.pi / 6, 1.0
    reps = 150
    means = _twin_means(l, c, n, th, reps=reps - 1, seed=79)
    dl = _delta(means)
    dirv = np.array([np.cos(th), np.sin(th)])
    lv = np.asarray(l)
    a_gamma = np.sqrt(lv) / (1.0 + lv) * g * dirv
    # trim keeping k=1 (score-regression form), same fixed-Q convention
    from simulator import spectrum

    rng = np.random.default_rng(81)
    p = int(round(c * n))
    r = 2
    Q, _ = np.linalg.qr(np.random.default_rng(79).standard_normal((p, r)))
    sv = np.sqrt(lv)
    acc = {g_: np.zeros(r) for g_ in (0.0, 1.0)}
    for g_ in (0.0, 1.0):
        for _ in range(reps):
            Xc, Yc = _draw_cell(Q, sv, r, n, th, rng, g_)
            d, V = spectrum(Xc)
            b = V.T @ (Xc.T @ Yc / n)
            mm = np.zeros(len(d))
            mm[:1] = 1.0 / d[:1]
            acc[g_] += Q.T @ (V @ (mm * b))
    dt = (acc[1.0] - acc[0.0]) / reps
    do = dl["ols"]
    # retained supercritical direction passes the F3 ideal scaled by the
    # finite-n BGN overlap (sample top axis captures xi of q1); dropped
    # subcritical direction passes ~0
    from de_formulas import bgn_overlap

    xi1 = bgn_overlap(l[0], c)
    pred_trim = np.array([xi1 * a_gamma[0], 0.0])
    assert np.linalg.norm(dt - pred_trim) < 0.05, (dt, pred_trim)
    assert abs(dt[1]) < 0.03, dt
    # OLS passes cap-law values on both coordinates
    cap = (1.0 + lv) / (c + lv)
    assert np.linalg.norm(do - cap * a_gamma) < 0.06, (do, cap * a_gamma)


def test_FROZEN_csv_anchor_ols_transmission():
    df = pd.read_csv(os.path.join(ROOT, "results",
                                  "estimation_cell_detail.csv"))
    row1 = df[(df.c == 0.2) & (df.profile == "sub") & (df.n == 500)
              & (np.isclose(df.theta, 0.5236)) & (df.r == 1)]
    row2 = df[(df.c == 5.0) & (df.profile == "sub") & (df.n == 2000)
              & (np.isclose(df.theta, 0.5236)) & (df.r == 1)]
    g = 1.0
    pred1 = np.sqrt(0.5 * np.sqrt(0.2)) / (1 + 0.5 * np.sqrt(0.2)) \
        * g * np.cos(0.5236)
    ll = 0.5 * np.sqrt(5.0)
    pred2 = np.sqrt(ll) / (5.0 + ll) * g * np.cos(0.5236)
    assert abs(pred1 - float(row1.ols_conf.iloc[0])) < 0.02
    assert abs(pred2 - float(row2.ols_conf.iloc[0])) < 0.02


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
