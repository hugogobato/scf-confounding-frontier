"""SCF Phase 2 runners (WP 2.1-2.4). Modes:
  correctness : pilot-schema estimation cells (ols + fixed-lambda ridges +
                pca_onatski) used for DE overlays
  estimation  : full frozen estimator roster (WP 2.2)
  crossover   : reduced roster over the g-grid strip (WP 2.2 figures)
  robustness  : reduced roster under DGP variants (WP 2.4)
  nullcal     : detection statistics only (WP 2.1 subset / WP 2.3 nulls)
  power       : detection statistics only, alternative cells (WP 2.3)
  alignment   : theta-sweep stats-only (WP 2.3 action 4)
  m2          : scalar-tau estimators (WP 2.4 M2 block)

Checkpoint/resume per cell exactly like the Phase 1 pilot (one parquet +
one means npz per config; finished cells skipped). All randomness flows from
simulator._rng so cells are reproducible and shardable in any grouping.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import json
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

from de_formulas import (
    bai_ng_select,
    eff_detect_spike,
    ledger_hash,
    onatski_select,
    ucm_strength,
)
from estimators import (
    RIDGE_LAM_GRID,
    _ridge_on_spectrum,
    eb_oracle_tau_factory,
    est_cevid_default,
    est_eb_cv_tau,
    est_eb_spectral,
    est_lava_transform_ols,
    est_pca_k,
    est_ridge_cv,
    fit_min_norm,
    fit_sdboost_linear_eb,
    make_oracle_estimator,
    center_columns,
    standardize_response,
)
from detection import compute_stats, probe_features, rejections
from simulator import Config, gen_data, spectrum

DOCS = Path(__file__).resolve().parents[1] / "docs"

ESTIMATION_FULL = [
    "ols", "ridge_cv", "pca_onatski", "pca_baing", "pca_oracle_r",
    "cevid_default", "lava_default", "sdboost_linear_eb", "eb_spectral",
    "eb_cv_tau", "eb_oracle_tau", "oracle_gamma",
]
ESTIMATION_REDUCED = [
    "ols", "ridge_cv", "pca_onatski", "pca_baing", "cevid_default",
    "sdboost_linear_eb", "eb_spectral", "eb_cv_tau", "eb_oracle_tau",
]
ROSTER_BY_MODE = {
    "estimation": ESTIMATION_FULL,
    "crossover": ESTIMATION_REDUCED,
    "robustness": ESTIMATION_REDUCED,
}


# ---------------------------------------------------------------------------
# rep-level runners
# ---------------------------------------------------------------------------


def _prepare(cfg: Config, rep: int):
    data = gen_data(cfg, rep)
    Xc = center_columns(data["X"])
    # RAW centered response for estimation metrics (pilot convention);
    # detection statistics standardize internally (scale-invariant).
    Yc = data["Y"] - data["Y"].mean()
    eig = spectrum(Xc)
    return data, Xc, Yc, eig


def _bias_dir(data) -> np.ndarray:
    b = data["Lam"] @ data["gam"]
    nb = float(np.linalg.norm(b))
    return b / nb if nb > 1e-12 else b


def _est_rows(cfg, Xc, Ys, eig, data, suffix: str, roster, rep: int):
    """Fit the frozen roster; returns list[dict] rows and {tag: diff vector}."""
    n, p = Xc.shape
    d, V = eig
    beta = data["beta"]
    u_b = _bias_dir(data)
    l_arr = np.asarray(cfg.l)
    gam = data["gam"]
    Ubasis = data["Lam"] / (cfg.sigma_u * np.sqrt(l_arr))[None, :]

    def pca_onatski(Xc_, Ys_, eig_, rng_=None):
        return est_pca_k(Xc_, Ys_, eig_, None, int(onatski_select(eig_[0]))), {}

    def pca_baing(Xc_, Ys_, eig_, rng_=None):
        nn, pp = Xc_.shape
        return est_pca_k(Xc_, Ys_, eig_, None,
                         int(bai_ng_select(eig_[0], nn, pp))), {}

    def pca_oracle(Xc_, Ys_, eig_, rng_=None):
        return est_pca_k(Xc_, Ys_, eig_, None, int(cfg.r)), {}

    rows, diffs = [], {}
    kmax_hint = None
    if getattr(cfg, "r_misspec", 0) != 0:
        kmax_hint = max(1, min(10, cfg.r + cfg.r_misspec))

    def eb_spectral_wrap(Xc_, Ys_, eig_, rng_=None):
        return est_eb_spectral(Xc_, Ys_, eig_, rng_, kmax_hint=kmax_hint)

    table = {
        "ols": lambda Xc_, Ys_, eig_, rng_: (fit_min_norm(Xc_, Ys_, eig_), {}),
        "ridge_cv": lambda Xc_, Ys_, eig_, rng_: (
            est_ridge_cv(Xc_, Ys_, eig_, rng_), {}),
        "pca_onatski": pca_onatski,
        "pca_baing": pca_baing,
        "pca_oracle_r": pca_oracle,
        "cevid_default": lambda Xc_, Ys_, eig_, rng_: (
            est_cevid_default(Xc_, Ys_, eig_), {}),
        "lava_default": lambda Xc_, Ys_, eig_, rng_: (
            est_lava_transform_ols(Xc_, Ys_, eig_), {}),
        "sdboost_linear_eb": fit_sdboost_linear_eb,
        "eb_spectral": eb_spectral_wrap,
        "eb_cv_tau": est_eb_cv_tau,
        "eb_oracle_tau": lambda Xc_, Ys_, eig_, rng_: (
            eb_oracle_tau_factory(l_arr, gam)(Xc_, Ys_, eig_), {}),
        "oracle_gamma": lambda Xc_, Ys_, eig_, rng_: (
            make_oracle_estimator(Ubasis, gam)(Xc_, Ys_, eig_), {}),
    }

    rows, diffs = [], {}
    for tag in roster:
        t0 = time.perf_counter()
        info = {}
        if tag.startswith("ridge_fixed|"):
            lam = float(tag.split("|")[1])
            bh = _ridge_on_spectrum(d, V, Xc, Ys, lam)
            info = {"lam": lam}
        else:
            bh, info = table[tag](Xc, Ys, eig, None)
        diff = bh - beta
        rows.append({
            "config_id": cfg.cid,
            "rep": rep,
            "seed": 20260823,
            "estimator": tag + suffix,
            "lambda": float(info.get("lam", info.get("tau", np.nan))),
            "k_select": int(info.get("m", -1)),
            "rel_err": float(np.linalg.norm(diff)),
            "proj_bias_dir": float(u_b @ diff),
            "tuning_grid": int(info.get("grid", -1)),
            "runtime_s": time.perf_counter() - t0,
        })
        diffs[tag + suffix] = diff.astype(np.float64)
        if tag == "ols" and suffix == "":
            det = compute_stats(Xc, Ys, eig, cfg.sigma_u)
            rows[-1].update({k: v for k, v in det.items() if k != "ktop"})
            rows[-1].update(rejections(det))
            rows[-1]["ucm_proxy"] = ucm_strength(d, n, p)
            rows[-1]["seff"] = eff_detect_spike(
                l_arr, gam, cfg.c, cfg.sigma_u ** 2)
    return rows, diffs


def run_estimation_rep(cfg: Config, rep: int, mode: str):
    roster = ROSTER_BY_MODE.get(mode, ESTIMATION_FULL)
    data, Xc, Ys, eig = _prepare(cfg, rep)
    rows, acc = _est_rows(cfg, Xc, Ys, eig, data, "", roster, rep)
    if cfg.twin_gamma0:
        Y0 = standardize_response(data["X"] @ data["beta"] + data["eps"])
        g0_data = {**data, "gam": np.zeros_like(data["gam"])}
        rows0, acc0 = _est_rows(cfg, Xc, Y0, eig, g0_data, "_g0", roster, rep)
        rows.extend(rows0)
        acc.update(acc0)
    return rows, acc


def run_correctness_rep(cfg: Config, rep: int):
    """Pilot-compatible schema: ols, ridge_fixed grid, pca_onatski."""
    roster = ["ols"] + [f"ridge_fixed|{lam}" for lam in RIDGE_LAM_GRID] \
        + ["pca_onatski"]
    data, Xc, Ys, eig = _prepare(cfg, rep)
    rows, acc = _est_rows(cfg, Xc, Ys, eig, data, "", roster, rep)
    if cfg.twin_gamma0:
        Y0 = standardize_response(data["X"] @ data["beta"] + data["eps"])
        g0_data = {**data, "gam": np.zeros_like(data["gam"])}
        rows0, acc0 = _est_rows(cfg, Xc, Y0, eig, g0_data, "_g0", roster, rep)
        rows.extend(rows0)
        acc.update(acc0)
    return rows, acc


def run_stats_rep(cfg: Config, rep: int):
    """Detection statistics only (nullcal / power / alignment modes)."""
    data, Xc, Yc, eig = _prepare(cfg, rep)
    det = compute_stats(Xc, Yc, eig, cfg.sigma_u)
    feats = probe_features(eig, det)
    row = {
        "config_id": cfg.cid, "rep": rep, "seed": 20260823,
        **det, **rejections(det),
        "ucm_proxy": ucm_strength(eig[0], cfg.n, cfg.p),
    }
    for i, v in enumerate(feats):
        row[f"f{i}"] = float(v)
    return [row], {}


def run_m2_rep(cfg: Config, rep: int):
    """Scalar-tau estimators under M2 (WP 2.4 treatment block)."""
    data, Xc, Yc, eig = _prepare(cfg, rep)
    d, V = eig
    n = cfg.n
    D = data["D"] - data["D"].mean()
    Ys = Yc
    tau_true = cfg.m2_tau

    def tau_from_joint(a_design):
        coef, *_ = np.linalg.lstsq(a_design, Ys, rcond=None)
        return float(coef[0])

    out = {}
    out["tau_ols"] = tau_from_joint(np.column_stack([D, Xc]))
    k = max(int(onatski_select(d)), cfg.r)
    S = Xc @ V[:, :k]
    out["tau_trim_onatski"] = tau_from_joint(np.column_stack([D, S]))
    S = Xc @ V[:, :cfg.r]
    out["tau_trim_oracle_r"] = tau_from_joint(np.column_stack([D, S]))
    lam = 1.0
    Sc = Xc.T @ Xc / n + lam * np.eye(len(d))
    rhs = Xc.T @ Ys / n
    c_vec = Xc.T @ D / n
    m = float(D @ D) / n
    Scinv_c = np.linalg.solve(Sc, c_vec)
    Scinv_rhs = np.linalg.solve(Sc, rhs)
    a_hat = (float(D @ Ys) / n - c_vec @ Scinv_rhs) / (m - c_vec @ Scinv_c)
    out["tau_ridge1"] = a_hat
    rows = []
    for tag, val in out.items():
        rows.append({
            "config_id": cfg.cid, "rep": rep, "seed": 20260823,
            "estimator": tag, "lambda": np.nan, "k_select": -1,
            "rel_err": abs(val - tau_true), "proj_bias_dir": np.nan,
            "tuning_grid": -1, "runtime_s": 0.0,
        })
    return rows, {}


MODE_RUNNER = {
    "correctness": run_correctness_rep,
    "estimation": lambda cfg, rep: run_estimation_rep(cfg, rep, "estimation"),
    "crossover": lambda cfg, rep: run_estimation_rep(cfg, rep, "crossover"),
    "robustness": lambda cfg, rep: run_estimation_rep(cfg, rep, "robustness"),
    "nullcal": run_stats_rep,
    "power": run_stats_rep,
    "alignment": run_stats_rep,
    "m2": run_m2_rep,
}


# ---------------------------------------------------------------------------
# cell runner with checkpoint/resume
# ---------------------------------------------------------------------------


def run_cell(job: dict):
    """job keys: config(dict), mode, reps, raw_path, means_path."""
    cfg = Config(**job["config"])
    mode = job["mode"]
    runner = MODE_RUNNER[mode]
    raw_path, means_path = Path(job["raw_path"]), Path(job["means_path"])
    if raw_path.exists():
        try:
            have = pd.read_parquet(raw_path, columns=["rep"])
            if have["rep"].nunique() >= job["reps"]:
                print(f"[skip] {mode} {cfg.profile}/{cfg.label} "
                      f"c={cfg.c:.2f} n={cfg.n} done", flush=True)
                return cfg.cid, 0.0
        except Exception:
            pass
    t0 = time.perf_counter()
    frames, acc = [], {}
    for rep in range(job["reps"]):
        rows, mean_diffs = runner(cfg, rep)
        frames.append(pd.DataFrame(rows))
        for k, v in mean_diffs.items():
            acc[k] = acc.get(k, 0.0) + v
    df = pd.concat(frames, ignore_index=True)
    df.to_parquet(raw_path, index=False)
    if acc:
        np.savez_compressed(
            means_path,
            **{k: (v / job["reps"]).astype(np.float32) for k, v in acc.items()},
        )
    dt = time.perf_counter() - t0
    print(f"[done] {mode} {cfg.profile}/{cfg.label} c={cfg.c:.2f} n={cfg.n} "
          f"theta={cfg.theta:.2f} -> {len(df)} rows in {dt:.1f}s", flush=True)
    return cfg.cid, dt


def run_jobs(jobs: list[dict], out_root, workers: int | None = None):
    out_root = Path(out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    for j in jobs:
        Path(j["raw_path"]).parent.mkdir(parents=True, exist_ok=True)
        Path(j["means_path"]).parent.mkdir(parents=True, exist_ok=True)
    if workers is None:
        workers = min(6, max(1, (os.cpu_count() or 4) - 2))
    results = [run_cell(j) for j in jobs] if workers <= 1 else \
        list(Pool(workers).imap_unordered(run_cell, jobs))
    meta = {
        "workers": workers, "jobs": len(jobs), "ledger_hash": ledger_hash(DOCS),
        "wall_s": sum(dt for _, dt in results),
        "per_cell_s": {cid: dt for cid, dt in results if dt > 0},
    }
    with open(out_root / "run_meta.json", "w") as f:
        json.dump(meta, f, indent=1)
    return meta
