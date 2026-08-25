"""T3(b) visibility-boundary falsifiers (docs/theory_T3_visibility_boundary.md).

Closed forms tested against fresh simulation of the exact probe objects:
  v0 = 1 + sigma_eps^2
  M0 = 1 + c(1 + sigma_eps^2)
  m0 = M0/v0
  E[q|g] = (M0 + g^2(omega + c))/(v0 + g^2),  omega = sum_j l_j dir_j^2
  saturation ceiling omega + c; kappa = omega_star - omega.
"""

import math
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "code"))
from simulator import Config, gen_data  # noqa: E402


def q_moments(c, prof, g, n_reps, n=1200, seed0=100):
    p = int(round(c * n))
    cfg = Config(n=n, p=p, r=3, l=(0.5 * math.sqrt(c),) * 3 if prof == "sub"
                 else (3 * math.sqrt(c), 0.5 * math.sqrt(c), 0.5 * math.sqrt(c)),
                 theta=math.pi / 6, g=g, profile=prof, label="t3b", q_fixed=True)
    qs = []
    for rep in range(n_reps):
        data = gen_data(cfg, rep)
        Xc = data["X"] - data["X"].mean(0, keepdims=True)
        Yc = data["Y"] - data["Y"].mean()
        s2 = float(Yc @ Yc) / len(Yc)
        b = Xc.T @ (Yc / math.sqrt(s2)) / len(Yc)
        qs.append(float(b @ b))
    return float(np.mean(qs))


def test_floor_and_shift_curve():
    c, prof, n_reps = 0.5, "sub", 60
    l_sub = 0.5 * math.sqrt(c)
    omega = l_sub * (math.cos(math.pi / 6) ** 2 + math.sin(math.pi / 6) ** 2)
    v0, M0 = 2.0, 1.0 + c * 2.0
    m0 = M0 / v0
    q0 = q_moments(c, prof, 0.0, n_reps)
    assert abs(q0 - m0) < 0.06 * m0, (q0, m0)
    g = 2.0
    q1 = q_moments(c, prof, g, n_reps)
    pred = (M0 + g * g * (omega + c)) / (v0 + g * g)
    assert abs(q1 - pred) < 0.05 * m0, (q1, pred)


def test_saturation_ceiling_and_kappa():
    """E[q|large g] approaches omega + c from the kappa side."""
    c, prof = 0.5, "sub"
    l_sub = 0.5 * math.sqrt(c)
    omega = l_sub
    ceiling = omega + c
    q_big = q_moments(c, prof, 6.0, 50)
    # still approaching: within 25% of the remaining gap after finite g
    m0 = (1.0 + 2.0 * c) / 2.0
    gap0 = m0 - ceiling
    assert abs(q_big - ceiling) < 0.75 * gap0
    kappa_pred = 0.5 - omega
    assert kappa_pred > 0 and abs((m0 - q_big) - kappa_pred) < 0.4 * kappa_pred


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
