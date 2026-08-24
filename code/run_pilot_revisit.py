"""WP 1.5 revisit pass: high-rep, twin-arm runs on the predeclared
decision-relevant cells + all auxiliary cells (research plan WP 1.5 trap
recovery: "1000 reps on the 6 most decision-relevant cells before the memo").

Fixes two issues discovered in the first pass:
  1. ridge mean-diff key collision (all lambdas overwrote one key) - fixed in
     simulator.py; revisit regenerates ridge overlays cleanly.
  2. the tiny-bias aux cells are Monte-Carlo-noise-limited; common-random-
     number twins (gamma vs gamma=0 on identical seeds) reduce the noise of
     the confounding-attributed bias estimate by orders of magnitude.

Outputs: data/pilot/revisit/*.parquet, data/pilot/revisit_means/*.npz,
data/pilot/revisit_summary.csv.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import math
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

from de_formulas import (
    eff_detect_spike,
    ledger_hash,
    minnorm_capture,
    ols_bias_vector,
    ridge_bias_vector,
    ridge_capture,
)
from run_pilot import LAM_GRID, N, R, profile_l
from simulator import Config, gamma_vector, gen_data, run_rep

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "pilot" / "revisit"
MEANS = ROOT / "data" / "pilot" / "revisit_means"
for d in (RAW, MEANS):
    d.mkdir(parents=True, exist_ok=True)

REPS = 1000


def build_grid() -> list[Config]:
    cfgs = []
    # six decision-relevant main cells at high reps (twins everywhere: paired
    # differences give clean gamma-attributed bias at any aspect ratio)
    picks = [
        ("sub", 0.2, 0.0), ("mixed", 0.2, 0.0), ("super", 0.2, 0.0),
        ("mixed", 0.8, 0.0), ("mixed", 5.0, 0.0), ("sub", 5.0, 0.0),
    ]
    for prof, c, th in picks:
        cfgs.append(Config(n=N, p=int(round(c * N)), r=R,
                           l=profile_l(prof, c), theta=th, profile=prof,
                           label="main", twin_gamma0=True, q_fixed=True,
                           reps=REPS))
    for c in (0.2, 0.8, 1.0, 5.0):
        for th in (math.pi / 6, math.pi / 2):
            cfgs.append(Config(n=N, p=int(round(c * N)), r=R,
                               l=profile_l("onespike", c), theta=th,
                               profile="onespike", label="aux",
                               twin_gamma0=True, q_fixed=True,
                               reps=400))
    return cfgs


def run_cell(cfg: Config):
    raw_path = RAW / f"{cfg.cid}.parquet"
    means_path = MEANS / f"{cfg.cid}.npz"
    if raw_path.exists() and means_path.exists():
        print(f"[skip] {cfg.profile}/{cfg.label} th={cfg.theta:.2f} c={cfg.c}",
              flush=True)
        return cfg.cid, 0.0
    t0 = time.perf_counter()
    acc: dict[str, np.ndarray] = {}
    frames = []
    for rep in range(cfg.reps):
        rows, mean_diffs = run_rep(cfg, rep, LAM_GRID)
        frames.append(pd.DataFrame(rows))
        for k, v in mean_diffs.items():
            acc[k] = acc.get(k, 0.0) + v
    df = pd.concat(frames, ignore_index=True)
    df.to_parquet(raw_path, index=False)
    np.savez_compressed(means_path,
                        **{k: (v / cfg.reps).astype(np.float32)
                           for k, v in acc.items()})
    print(f"[done] {cfg.profile}/{cfg.label} th={cfg.theta:.2f} c={cfg.c} "
          f"-> {len(df)} rows in {(time.perf_counter()-t0)/60:.1f} min",
          flush=True)
    return cfg.cid, time.perf_counter() - t0


def summarize(cfgs: list[Config]) -> pd.DataFrame:
    recs = []
    for cfg in cfgs:
        mp = MEANS / f"{cfg.cid}.npz"
        rp = RAW / f"{cfg.cid}.parquet"
        if not mp.exists():
            continue
        means = np.load(mp)
        df = pd.read_parquet(rp) if rp.exists() else None
        gam = gamma_vector(cfg)
        l_arr = np.asarray(cfg.l)

        sim_total = float(np.linalg.norm(means["ols"].astype(float)))
        sim_conf = float(np.linalg.norm(
            (means["ols"] - means["ols_g0"]).astype(float)))
        if cfg.c <= 1.0:
            coef = ols_bias_vector(l_arr, gam)
        else:
            cap = minnorm_capture(l_arr, cfg.c)
            coef = cap * ols_bias_vector(l_arr, gam)
        pred_vec = gen_data(cfg, 0)["Q"] @ coef
        pred_conf = float(np.linalg.norm(pred_vec))
        pred_total = pred_conf
        # directional estimate: the ||mean-diff|| norm carries a
        # sqrt(p/reps) Monte-Carlo noise floor (measured ~0.025 at
        # c = 0.2, 400 reps) that buries 0.01-scale signals; the
        # projection on the predicted direction is the sharp estimator
        md = (means["ols"] - means["ols_g0"]).astype(float)
        proj_pred = float(md @ (pred_vec / pred_conf))

        ridge_rows = []
        for lam in LAM_GRID:
            key = f"ridge_fixed@{lam}"
            gkey = f"ridge_fixed_g0@{lam}"
            if key not in means:
                continue
            sim_r = float(np.linalg.norm(
                (means[key] - (means[gkey] if gkey in means else 0))
                .astype(float)))
            if cfg.c <= 1.0:
                pred_r = float(np.linalg.norm(
                    ridge_bias_vector(l_arr, gam, lam)))
            else:
                capr = ridge_capture(l_arr, lam, cfg.c)
                coef = capr * np.sqrt(l_arr * cfg.sigma_u ** 2) / \
                    (cfg.sigma_u ** 2 * (1.0 + l_arr) + lam) * gam
                pred_r = float(np.linalg.norm(coef))
            ridge_rows.append({"lambda": lam, "sim": sim_r, "pred": pred_r})

        ols_df = df[df["estimator"] == "ols"] if df is not None else None
        recs.append({
            "label": cfg.label, "profile": cfg.profile,
            "theta_deg": round(math.degrees(cfg.theta), 1), "c": cfg.c,
            "reps": cfg.reps,
            "bias_ratio_ols_sim": sim_total,
            "bias_ratio_conf_sim_twin": sim_conf,
            "bias_ratio_pred": pred_total,
            "proj_pred_sim": proj_pred,
            "dev_total_pct": 100 * abs(sim_total - pred_total) / max(pred_total, 1e-12),
            "dev_conf_pct": 100 * abs(sim_conf - pred_conf) / max(pred_conf, 1e-12),
            "dev_proj_pct": 100 * abs(proj_pred - pred_conf) / max(pred_conf, 1e-12),
            "outlier_rate_tw99": float(ols_df["outlier99"].mean())
            if ols_df is not None else np.nan,
            "lam_max_sim_mean": float(ols_df["lam_max_cov"].mean())
            if ols_df is not None else np.nan,
            "xi1_sim": float(ols_df["overlap1"].dropna().mean())
            if ols_df is not None else np.nan,
            "ridge_curves": ridge_rows,
            "seff": eff_detect_spike(l_arr, gam, cfg.c, cfg.sigma_u),
            "ledger_hash": ledger_hash(ROOT / "docs"),
        })
    return pd.DataFrame(recs)


if __name__ == "__main__":
    t0 = time.perf_counter()
    cfgs = build_grid()
    print(f"revisit grid: {len(cfgs)} cells x {REPS} reps", flush=True)
    workers = min(6, max(1, (os.cpu_count() or 4) - 2))
    with Pool(workers) as pool:
        pool.map(run_cell, cfgs)
    s = summarize(cfgs)
    s.to_csv(ROOT / "data" / "pilot" / "revisit_summary.csv", index=False)
    parts = [pd.read_parquet(RAW / f"{c.cid}.parquet") for c in cfgs
             if (RAW / f"{c.cid}.parquet").exists()]
    pd.concat(parts, ignore_index=True).to_parquet(
        ROOT / "data" / "pilot" / "revisit_results.parquet", index=False)
    print(f"revisit done in {(time.perf_counter()-t0)/60:.1f} min wall",
          flush=True)
