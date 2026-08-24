"""Phase 2 cost-calibration micro-pilot (plan Section 10.1 rule 1).

Measures per-rep wall time for each mode at representative (n, p) points on
the CURRENT machine, writes data/sim/cost_model.json used by make_shards.py
to group cells into <= 5 h shards. Small reps only; do not run long.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import json
import math
import time
from pathlib import Path

import numpy as np

from simulator import Config
from runners import (
    run_correctness_rep,
    run_estimation_rep,
    run_stats_rep,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "sim"
OUT.mkdir(parents=True, exist_ok=True)

POINTS = [
    # (n, p, r, prof) covering the grid's cost corners
    (500, 50, 1, "sub"), (500, 5000, 25, "mixed"),
    (2000, 200, 1, "sub"), (2000, 2000, 25, "mixed"),
    (2000, 10000, 5, "mixed"), (2000, 20000, 5, "mixed"),
    (8000, 800, 1, "sub"), (8000, 16000, 5, "mixed"),
]
REPS_TIMED = 3


def time_mode(fn, cfg, reps):
    ts = []
    for rep in range(reps):
        t0 = time.perf_counter()
        fn(cfg, rep)
        ts.append(time.perf_counter() - t0)
    return float(np.median(ts)), float(np.max(ts))


def main():
    model = {"machine": "local i9-13900H under ambient load", "points": []}
    for n, p, r, prof in POINTS:
        sc = math.sqrt(p / n)
        l = tuple([3.0 * sc] + [0.5 * sc] * (r - 1)) if prof == "mixed" \
            else tuple([0.5 * sc] * r)
        kw = dict(n=n, p=p, r=r, l=l, theta=math.pi / 6, g=1.0,
                  profile=prof, label="cost", q_fixed=True)
        est_cfg = Config(**kw)
        stats_cfg = Config(**{**kw, "g": 0.0})
        try:
            med_est, max_est = time_mode(run_correctness_rep, est_cfg, REPS_TIMED)
            med_full, _ = time_mode(
                lambda c_, r_: run_estimation_rep(c_, r_, "estimation"),
                est_cfg, 2) if n <= 2000 else (float("nan"), float("nan"))
            med_stats, max_stats = time_mode(run_stats_rep, stats_cfg, REPS_TIMED)
        except MemoryError:
            print(f"skip (n={n}, p={p}): MemoryError")
            continue
        rec = {
            "n": n, "p": p, "r": r, "profile": prof,
            "s_per_rep_correctness": round(med_est, 3),
            "s_per_rep_estimation_full": round(med_full, 3),
            "s_per_rep_statsonly": round(med_stats, 3),
        }
        model["points"].append(rec)
        print(rec, flush=True)
    with open(OUT / "cost_model.json", "w") as f:
        json.dump(model, f, indent=1)
    print("->", OUT / "cost_model.json")


if __name__ == "__main__":
    main()
