"""Deterministic smoke run (WP 1.4 action 3): end-to-end in under 1 minute,
writing a valid parquet at data/smoke/run_smoke.parquet."""
from __future__ import annotations

import math
import os
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))

from simulator import Config, run_rep


def test_smoke_end_to_end(tmp_path):
    t0 = time.perf_counter()
    cfg = Config(n=500, p=200, r=2, l=(1.5, 0.5), theta=math.pi / 6,
                 twin_gamma0=True, profile="smoke", label="smoke")
    rows, _ = run_rep(cfg, 0, (0.1, 1.0, 10.0))
    df = pd.DataFrame(rows)
    out = tmp_path / "run_smoke.parquet"
    df.to_parquet(out, index=False)
    back = pd.read_parquet(out)
    assert len(back) == len(df) and len(back) > 0
    assert {"config_id", "rep", "estimator", "rel_err"}.issubset(back.columns)
    assert time.perf_counter() - t0 < 60.0
