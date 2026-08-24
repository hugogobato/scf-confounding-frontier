"""WP 1.5 pilot runner: builds the predeclared grid and runs it in parallel.

Grid (research plan WP 1.5 + predeclared amendments documented in
docs/pilot_memo.md):
  main : c in {0.2, 1.0, 5.0} x profile in {sub, mixed, super} x theta in
         {0, pi/6, pi/2}, r=3, n=2000, p=round(c*n), g=1, 200 reps.
  aux  : profile "onespike" (l = 3*sqrt(c), 0.01, 0.01) x theta in
         {pi/6, pi/2} x same c's (visible-spike / negligible-bias boundary).
Twin gamma=0 arms run for c > 1 cells (where min-norm shrinkage contaminates
total bias); at c <= 1 the exact identity makes confounding bias equal total
bias, so twins would be redundant.

Checkpointing: one parquet + one means-npz per cell; finished cells are
skipped on restart (plan Section 10.1 rule 4).
"""
from __future__ import annotations

import os

# Cap BLAS threads BEFORE numpy is imported anywhere in this process tree
# (research plan Section 10.1: worker-level single-thread BLAS, <= 8 workers).
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import json
import math
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

from de_formulas import (
    eff_detect_spike,
    ledger_hash,
    minnorm_bias_vector,
    minnorm_total_bias_norm,
    ols_bias_vector,
    ridge_bias_vector,
    ridge_capture,
)
from simulator import Config, gamma_vector, run_rep

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "pilot" / "raw"
MEANS = ROOT / "data" / "pilot" / "means"
DOCS = ROOT / "docs"
for d in (RAW, MEANS):
    d.mkdir(parents=True, exist_ok=True)

N = 2000
REPS = 200
R = 3
LAM_GRID = (0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0)


def profile_l(name: str, c: float) -> tuple[float, ...]:
    sc = math.sqrt(c)
    if name == "sub":
        return (0.5 * sc,) * R
    if name == "mixed":
        return (3.0 * sc, 0.5 * sc, 0.5 * sc)
    if name == "super":
        return (3.0 * sc,) * R
    if name == "onespike":
        # secondary spikes negligible enough that sqrt(l)/(1+l) <= 0.01,
        # i.e., bias ratio <= 0.02 with g = 1 (visible-outlier / harmless cell)
        return (3.0 * sc, 1e-4, 1e-4)
    raise ValueError(name)


def build_grid() -> list[Config]:
    # q_fixed=True: loading directions are drawn once per config so that
    # mean(beta_hat) - beta estimates the bias CONDITIONAL on the loading
    # geometry (the object the DE formulas predict). With Haar-random Q per
    # rep the mean bias vector would average to zero by symmetry.
    cfgs = []
    # c = 1.0 is kept for ridge/eigenstructure overlays but its OLS mean-bias
    # metric is UNSTABLE (E[(X'X)^{-1}] diverges at aspect 1; measured
    # mean-ratios of order 30-165 across reps). A c = 0.8 block replaces it
    # for the p < n bias-overlay claims (documented deviation, pilot memo).
    for c in (0.2, 0.8, 1.0, 5.0):
        p = int(round(c * N))
        for prof in ("sub", "mixed", "super"):
            for th in (0.0, math.pi / 6, math.pi / 2):
                cfgs.append(Config(
                    n=N, p=p, r=R, l=profile_l(prof, c), theta=th,
                    profile=prof, label="main",
                    twin_gamma0=(c > 1.0), q_fixed=True, reps=REPS))
        for th in (math.pi / 6, math.pi / 2):
            cfgs.append(Config(
                n=N, p=p, r=R, l=profile_l("onespike", c), theta=th,
                profile="onespike", label="aux",
                twin_gamma0=(c > 1.0), q_fixed=True, reps=REPS))
    return cfgs


def run_cell(cfg: Config):
    raw_path = RAW / f"{cfg.cid}.parquet"
    means_path = MEANS / f"{cfg.cid}.npz"
    expected_rows_per_rep = len(LAM_GRID) + 3  # ols + ridges + pca_onatski
    expected = cfg.reps * expected_rows_per_rep * (2 if cfg.twin_gamma0 else 1)
    if raw_path.exists():
        try:
            have = pd.read_parquet(raw_path, columns=["rep"])
            if len(have) >= expected:
                print(f"[skip] {cfg.profile}/{cfg.label} theta={cfg.theta:.2f} "
                      f"c={cfg.c} done ({len(have)} rows)", flush=True)
                return cfg.cid, 0.0
        except Exception:
            pass

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
    np.savez_compressed(
        means_path,
        **{k: (v / cfg.reps).astype(np.float32) for k, v in acc.items()},
    )
    dt = time.perf_counter() - t0
    print(f"[done] {cfg.profile}/{cfg.label} theta={cfg.theta:.2f} c={cfg.c} "
          f"-> {len(df)} rows in {dt:.1f}s", flush=True)
    return cfg.cid, dt


def summarize(cfgs: list[Config]) -> pd.DataFrame:
    recs = []
    for cfg in cfgs:
        cid = cfg.cid
        mp = MEANS / f"{cid}.npz"
        rp = RAW / f"{cid}.parquet"
        if not mp.exists() or not rp.exists():
            continue
        means = np.load(mp)
        df = pd.read_parquet(rp)
        ols_df = df[df["estimator"] == "ols"]
        gam = gamma_vector(cfg)
        l_arr = np.asarray(cfg.l)
        sup = l_arr > math.sqrt(cfg.c)

        # OLS bias ratios. For c > 1 the fit-artifact term (1/c - 1) beta_perp
        # cancels in mean(beta_hat - beta) because beta is redrawn per rep
        # (A4a, E[beta] = 0); hence the simulated total-bias ratio and the
        # gamma-attributed twin difference estimate the SAME functional and
        # the prediction is the confounding capture vector (F8).
        sim_total = float(np.linalg.norm(means["ols"].astype(float)))
        if cfg.c <= 1.0:
            pred_total = float(np.linalg.norm(ols_bias_vector(l_arr, gam)))
            pred_conf = pred_total
            sim_conf = sim_total
        else:
            conf_vec, _ = minnorm_bias_vector(
                l_arr, gam, cfg.c,
                np.zeros(cfg.r))  # beta coords irrelevant: fit term cancels
            pred_conf = float(np.linalg.norm(conf_vec))
            pred_total = pred_conf
            sim_conf = (float(np.linalg.norm(
                means["ols"].astype(float) - means["ols_g0"].astype(float)))
                if "ols_g0" in means else np.nan)

        # empirical capture coefficients (fixed-Q cells only): projection of
        # the mean bias vector on u_j divided by the population coefficient
        # times gamma_j (the bias along u_j is cap_j * coef_j * gamma_j)
        emp_cap = []
        if cfg.q_fixed:
            from simulator import gen_data
            Qref = gen_data(cfg, 0)["Q"]
            md = means["ols"].astype(float)
            for j in range(cfg.r):
                denom = math.sqrt(l_arr[j] / cfg.sigma_u ** 2) \
                    / (1.0 + l_arr[j]) * gam[j]
                emp_cap.append(float(Qref[:, j] @ md / denom)
                               if abs(denom) > 1e-12 else np.nan)

        # RMS-error functional (secondary validation target of F8): the
        # min-norm shrinkage artifact lives in per-rep ERROR, not in bias
        rms_err_sim = float(np.sqrt((ols_df["rel_err"] ** 2).mean()))
        pred_rms_err = (minnorm_total_bias_norm(l_arr, gam, cfg.c, cfg.p)
                        if cfg.c > 1.0 else np.nan)

        # ridge overlays
        ridge_devs = []
        for lam in LAM_GRID:
            key = f"ridge_fixed|{lam}"
            if key not in means:
                continue
            sim_r = float(np.linalg.norm(means[key].astype(float)))
            if cfg.c <= 1.0:
                pred_r = float(np.linalg.norm(ridge_bias_vector(l_arr, gam, lam)))
            else:
                capr = ridge_capture(l_arr, lam, cfg.c)
                coef = capr * np.sqrt(l_arr * cfg.sigma_u ** 2) / \
                    (cfg.sigma_u ** 2 * (1.0 + l_arr) + lam) * gam
                pred_r = float(np.linalg.norm(coef))
            ridge_devs.append(abs(sim_r - pred_r) / max(pred_r, 1e-12))
        # detection stats from raw
        outlier_rate = float(ols_df["outlier99"].mean())
        bbp_sim = float(ols_df["lam_max_cov"].mean())
        xi1_sim = (float(ols_df["overlap1"].dropna().mean())
                   if ols_df["overlap1"].notna().any() else np.nan)
        seff = eff_detect_spike(l_arr, gam, cfg.c, cfg.sigma_u)

        recs.append({
            "config_id": cid,
            "label": cfg.label,
            "profile": cfg.profile,
            "theta_deg": round(math.degrees(cfg.theta), 1),
            "c": cfg.c,
            "n": cfg.n,
            "p": cfg.p,
            "g": cfg.g,
            "l1": float(l_arr[0]),
            "l2": float(l_arr[1]),
            "l3": float(l_arr[2]),
            "twin_gamma0": cfg.twin_gamma0,
            "stable_mean_ols": bool(cfg.c < 1.0),
            "reps_done": int(len(ols_df)),
            "bias_ratio_ols_sim": sim_total,
            "bias_ratio_ols_pred": pred_total,
            "dev_ols_pct": 100.0 * abs(sim_total - pred_total) / max(pred_total, 1e-12),
            "bias_ratio_conf_sim": sim_conf,
            "bias_ratio_conf_pred": pred_conf,
            "emp_cap1": emp_cap[0] if len(emp_cap) > 0 else np.nan,
            "emp_cap2": emp_cap[1] if len(emp_cap) > 1 else np.nan,
            "emp_cap3": emp_cap[2] if len(emp_cap) > 2 else np.nan,
            "rms_err_sim": rms_err_sim,
            "pred_rms_err": pred_rms_err,
            "outlier_rate_tw99": outlier_rate,
            "lam_max_sim_mean": bbp_sim,
            "bbp_pred": float(ols_df["bbp_pred"].iloc[0]),
            "xi1_sim": xi1_sim if not np.isnan(xi1_sim) else np.nan,
            "xi1_pred": float(ols_df["xi1_pred"].iloc[0]),
            "seff": seff,
            "ridge_dev_mean_pct": 100.0 * float(np.nanmean(ridge_devs)),
            "ledger_hash": ledger_hash(DOCS),
        })
    return pd.DataFrame(recs)


def main():
    t0 = time.perf_counter()
    cfgs = build_grid()
    print(f"grid: {len(cfgs)} cells x {REPS} reps; ledger hash "
          f"{ledger_hash(DOCS)}", flush=True)
    # research plan Section 10.1: cap workers below physical core count and
    # leave headroom when other experiments are running (checked via uptime)
    workers = min(6, max(1, (os.cpu_count() or 4) - 2))
    with Pool(workers) as pool:
        for cid, dt in pool.imap_unordered(run_cell, cfgs):
            pass
    summary = summarize(cfgs)
    out_csv = ROOT / "data" / "pilot" / "cell_summary.csv"
    summary.to_csv(out_csv, index=False)
    # consolidated parquet (one row per config x rep x estimator)
    parts = [pd.read_parquet(RAW / f"{c.cid}.parquet") for c in cfgs
             if (RAW / f"{c.cid}.parquet").exists()]
    full = pd.concat(parts, ignore_index=True)
    full.to_parquet(ROOT / "data" / "pilot" / "pilot_results.parquet", index=False)
    print(f"all cells done in {(time.perf_counter() - t0) / 60:.1f} min wall; "
          f"{len(full)} rows -> pilot_results.parquet; summary -> {out_csv}",
          flush=True)
    with open(ROOT / "data" / "pilot" / "run_meta.json", "w") as f:
        json.dump({"workers": workers, "reps": REPS, "n_cells": len(cfgs),
                   "wall_min": (time.perf_counter() - t0) / 60,
                   "ledger_hash": ledger_hash(DOCS)}, f, indent=1)


if __name__ == "__main__":
    main()
