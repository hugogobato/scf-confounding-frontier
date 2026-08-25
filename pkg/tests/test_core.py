import numpy as np

from confounderalarm import fit_alarm


def _synthetic(n=320, p=140, r=2, l=(6.0, 6.0), g=0.0, seed=7):
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.standard_normal((p, r)))
    Lam = Q * np.sqrt(np.asarray(l))[None, :]
    f = rng.standard_normal((n, r))
    X = f @ Lam.T + rng.standard_normal((n, p))
    beta = rng.standard_normal(p)
    beta /= np.linalg.norm(beta)
    gam = g * np.ones(r) / np.sqrt(r)
    Y = X @ beta + f @ gam + rng.standard_normal(n)
    return Y, X


def test_no_alarm_on_clean_design():
    Y, X = _synthetic(g=0.0, seed=11)
    rep = fit_alarm(Y, X, n_perm=200, seed=3)
    assert rep.p_value > 0.01
    assert not rep.alarm


def test_alarm_fires_on_dense_confounding():
    Y, X = _synthetic(g=4.0, seed=12)
    rep = fit_alarm(Y, X, n_perm=200, seed=3)
    assert rep.alarm
    assert rep.p_value <= 0.05
    assert rep.g_star is not None and np.isfinite(rep.g_star)


def test_adjustment_recovers_tau():
    # c > 1 with separated spikes replicates the Phase-2 M2 regime where
    # raw-OLS tau error inflates (min-norm capture artifacts + confounding)
    # while Onatski trim-then-regress stays flat on tau_true = 1.
    rng = np.random.default_rng(5)
    n, p, r = 300, 600, 2
    Q, _ = np.linalg.qr(rng.standard_normal((p, r)))
    Lam = Q * np.sqrt(np.array([16.0, 4.0]))[None, :]
    f = rng.standard_normal((n, r))
    X = f @ Lam.T + rng.standard_normal((n, p))
    kpi = max(3, p // 100)
    pi = np.zeros(p)
    pi[:kpi] = 1.0 / np.sqrt(kpi)
    delta = rng.standard_normal(r)
    delta *= 0.3 / np.linalg.norm(delta)
    D = X @ pi + f @ delta + rng.standard_normal(n)
    beta = rng.standard_normal(p)
    beta /= np.linalg.norm(beta)
    gam = 2.0 * np.ones(r) / np.sqrt(r)
    Y = 1.0 * D + X @ beta + f @ gam + rng.standard_normal(n)
    rep = fit_alarm(Y, X, D=D, n_perm=150, seed=9)
    adj = rep.adjustment
    assert adj["method"] == "onatski_trim"
    assert abs(adj["tau_trim"] - 1.0) < abs(adj["tau_ols"] - 1.0)
    assert abs(adj["tau_trim"] - 1.0) < 0.15
    assert abs(np.sign(adj["tau_trim"] - 0.0) - 1.0) < 1e-12
