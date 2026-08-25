"""SCF Phase 3 benchmark machinery (WP 3.2/3.3): semi-synthetic injection on
real designs, calibrated alarm + mandatory baselines, arms runner with
checkpoint/resume, and the F12-law frontier prediction g* per benchmark.

Frozen spec: configs/benchmarks_frozen.yaml (v1). Two-pass order is binding:
pass 1 = null + perm_null arms only (thresholds and g* into
results/benchmark_freeze.json); pass 2 = everything else. All randomness
flows from SeedSequence([GLOBAL_SEED, cell_hash, rep]) so matched-null twins
share beta/f/eps draws with positive arms at equal rep index.

Injection model (identical semantics to simulator.gen_data M1/M2):
    X_obs = Xc_base + f @ Lam'
    Y     = X_obs @ beta + f @ gam + eps
so Cov(X_obs, Y) = Sigma_obs beta + Lambda gamma with
Sigma_obs = Sigma_base + Lambda Lambda', exactly the M1 algebra; the ground
truth is beta wrt the OBSERVED design. Null twins drop f from BOTH blocks.

Baselines implemented here (both flagged APPROXIMATE transcriptions, same
policy as Phase 2's ucm_strength):
  ucm_rho  : response-aware confounding-variance share proxy in the spirit
             of Rendsburg et al. (2022): sum_j g2hat_j l_j/(1+l_j) / mean(d)
             over the alarm coordinates, g2hat from estimators.estimate_gamma2.
  js_asym  : Janzing-Schoelkopf-style spectral asymmetry (2018): eigenvalue
             drops of the design covariance after removing the rank-one
             response-explained component, max relative drop over the top-K
             interlacing roots.
Both are calibrated by the permuted-Y null arm (frozen q95 thresholds), not
by their papers' asymptotics.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    import os

    os.environ.setdefault(_v, "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import pandas as pd  # noqa: E402

from de_formulas import (  # noqa: E402
    bbp_location,
    ledger_hash,
    minnorm_capture,
    mp_edges,
    onatski_select,
    tw_mu_sigma,
    tw_threshold,
)
from detection import compute_stats, rejections  # noqa: E402
from estimators import estimate_gamma2  # noqa: E402
from simulator import GLOBAL_SEED, spectrum  # noqa: E402

BENCH_DIR = ROOT / "data" / "benchmarks"
RAW_DIR = BENCH_DIR / "raw"
FREEZE_PATH = ROOT / "results" / "benchmark_freeze.json"


# ---------------------------------------------------------------------------
# benchmark objects
# ---------------------------------------------------------------------------


class Bench:
    """One processed real design plus its frozen spectral profile."""

    def __init__(self, name: str):
        self.name = name
        z = np.load(BENCH_DIR / f"{name}.npz", allow_pickle=True)
        self.X = np.ascontiguousarray(z["X"], dtype=np.float64)
        meta = json.loads(str(z["meta_json"]))
        self.meta = meta
        self.n, self.p = self.X.shape
        self.c = self.p / self.n
        Xc = self.X - self.X.mean(axis=0, keepdims=True)
        self.Xc = Xc
        self.eig = spectrum(Xc)
        d = self.eig[0]
        self.d = d
        self.se2 = float(max(np.quantile(d, 0.25),
                             1e-3 * float(np.mean(d))))
        self.ktop = int(min(max(onatski_select(d), 1), 10))
        self.r_inj = self.ktop
        self.l_hat = np.maximum(d[:10] / self.se2 - 1.0, 0.0)
        self.V = self.eig[1]

    def lam_matrix(self, r: int) -> np.ndarray:
        r = int(max(1, min(r, 10)))
        sds = np.sqrt(np.maximum(self.d[:r] - self.se2, 1e-12))
        return self.V[:, :r] * sds[None, :]


def cid_for(config_name: str, arm: str) -> str:
    import hashlib

    payload = f"benchmarks_frozen_v1|{config_name}|{arm}"
    return hashlib.sha256(payload.encode()).hexdigest()[:12]


def rng_for(config_name: str, arm: str, rep: int) -> np.random.Generator:
    return np.random.default_rng(
        np.random.SeedSequence(
            entropy=[GLOBAL_SEED, int(cid_for(config_name, arm)[:8], 16), rep]
        )
    )


def gamma_dir(bench: Bench, kind: str, r: int) -> np.ndarray:
    dirv = np.zeros(max(r, bench.r_inj))
    if kind == "spread":
        dirv[:r] = 1.0 / np.sqrt(max(r, 1))
    elif kind == "top":
        dirv[0] = 1.0
    elif kind == "weak":
        dirv[r - 1] = 1.0
    else:
        raise ValueError(kind)
    return dirv[:r]


# ---------------------------------------------------------------------------
# one replication
# ---------------------------------------------------------------------------


def _draw_core(bench: Bench, config_name: str, arm: str, rep: int,
               r: int):
    rng = rng_for(config_name, arm, rep)
    beta = rng.standard_normal(bench.p)
    beta /= float(np.linalg.norm(beta))
    f = rng.standard_normal((bench.n, r))
    eps = rng.standard_normal(bench.n)
    pi = None
    delta = None
    nu = None
    wants_m2 = bool(arm_specs()["configs"].get(config_name, {}).get("m2"))
    if wants_m2:
        kpi = max(3, bench.p // 100)
        pi = np.zeros(bench.p)
        pi[:kpi] = 1.0 / np.sqrt(kpi)
        delta = rng.standard_normal(r)
        delta *= 0.3 / max(float(np.linalg.norm(delta)), 1e-12)
        nu = rng.standard_normal(bench.n)
    return dict(beta=beta, f=f, eps=eps, pi=pi, delta=delta, nu=nu)


def run_rep(bench: Bench, config_name: str, arm_spec: dict, arm: str,
            rep: int, g_star: float | None,
            coord_scales: np.ndarray | None = None):
    """Return (rows, mean_diff_ols or None)."""
    t0 = time.perf_counter()
    kind = arm_spec["type"]
    r = bench.r_inj + int(arm_spec.get("r_delta", 0))
    r = int(min(max(r, 1), 10)) if kind != "split_half_null" else bench.r_inj
    core = _draw_core(bench, config_name, arm, rep, r)
    beta, f, eps = core["beta"], core["f"], core["eps"]
    hetero = bool(arm_spec.get("hetero_eps", False))
    if hetero:
        w = rng_for(config_name, arm + "_het", rep).chisquare(1.0, bench.n)
        eps = eps * np.sqrt(w / w.mean())
    gscale = float(arm_spec.get("g_scale", 1.0))
    gam_scale = 0.0 if arm_spec.get("gamma_zero") else (
        (g_star if g_star is not None else 1.0) * gscale)
    dirv = gamma_dir(bench, arm_spec.get("dir", "spread"), r)
    gam = gam_scale * dirv

    rows = []
    mean_diff = None

    def stats_block(Xc_used, Y_used, eig_used, tag_rows):
        det = compute_stats(Xc_used, Y_used - Y_used.mean(), eig_used)
        zeta = raw_z_coords(Xc_used, Y_used - Y_used.mean(), eig_used,
                            bench.ktop)
        t_alarm, k_used = bench_alarm_stat(Xc_used, Y_used - Y_used.mean(),
                                           eig_used, bench, coord_scales)
        row = {
            "config_id": cid_for(config_name, arm),
            "config_name": config_name, "arm": arm, "rep": rep,
            **det, **rejections(det),
            "T_bench": t_alarm,
            "ktop_bench": k_used,
            "ucm_rho": ucm_rho(Xc_used, Y_used - Y_used.mean(), eig_used,
                               bench),
            "js_asym": js_asym_from_stats(eig_used, Xc_used,
                                          Y_used - Y_used.mean()),
        }
        for j, zv in enumerate(zeta):
            row[f"zeta{j}"] = float(zv)
        tag_rows.append(row)
        return row

    if kind == "bench":
        Lam = bench.lam_matrix(r)
        Xobs = bench.Xc + f @ Lam.T
        Y = Xobs @ beta + f @ gam + eps
        eig_obs = spectrum(Xobs)
        Xc_obs = Xobs - Xobs.mean(axis=0, keepdims=True)
        row = stats_block(Xc_obs, Y, eig_obs, rows)
        b_ols = _ols_beta(Xc_obs, Y - Y.mean(), eig_obs)
        mean_diff = b_ols - beta
        row["rel_bias_dir"] = float(np.linalg.norm(mean_diff))
    elif kind == "permute_y_of_null":
        Xnull = bench.Xc
        Y0 = Xnull @ beta + eps
        prng = rng_for(config_name, arm + "_p", rep)
        Yperm = prng.permutation(Y0)
        stats_block(Xnull, Yperm, bench.eig, rows)
    elif kind == "split_half_null":
        prng = rng_for(config_name, arm + "_s", rep)
        idx = np.sort(prng.choice(bench.n, size=bench.n // 2, replace=False))
        Xh = bench.Xc[idx]
        Yh = Xh @ beta + eps[idx]
        Xh = Xh - Xh.mean(axis=0, keepdims=True)
        eig_h = spectrum(Xh)
        stats_block(Xh, Yh, eig_h, rows)
    elif kind == "m2":
        Lam = bench.lam_matrix(r)
        Xobs = bench.Xc + f @ Lam.T
        D = Xobs @ core["pi"] + f @ core["delta"] + core["nu"]
        Y = 1.0 * D + Xobs @ beta + f @ gam + eps
        eig_obs = spectrum(Xobs)
        Xc_obs = Xobs - Xobs.mean(axis=0, keepdims=True)
        Dc = D - D.mean()
        Yc = Y - Y.mean()
        taus = tau_estimators(Xc_obs, Dc, Yc, eig_obs, bench)
        row = stats_block(Xc_obs, Y, eig_obs, rows)
        row.update(taus)
    else:
        raise ValueError(kind)

    for row in rows:
        row["runtime_s"] = time.perf_counter() - t0
    return rows, mean_diff


def _ols_beta(Xc, Yc, eig):
    from simulator import fit_ols

    return fit_ols(Xc, Yc, eig)


def ucm_rho(Xc: np.ndarray, Yc: np.ndarray, eig, bench: Bench) -> float:
    """Response-aware confounding-variance share proxy (APPROXIMATE UCM).

    rho_hat = sum_{j<ktop} g2hat_j l_j/(1+l_j) / mean(d): the estimated
    response-linked variance carried by the leading directions, divided by
    the total design variance share of the bulk. Monotone in the injected
    link by construction; calibration is permutation-based upstream."""
    n, p = Xc.shape
    c = p / n
    k = min(bench.ktop, len(bench.l_hat))
    l_hat = bench.l_hat[:k]
    g2 = estimate_gamma2(Xc, Yc, eig, l_hat, c)
    conf_var = float(np.sum(g2 * l_hat / (1.0 + l_hat)))
    return float(conf_var / max(float(np.mean(eig[0])), 1e-12))


def raw_z_coords(Xc: np.ndarray, Yc: np.ndarray, eig,
                 ktop: int) -> np.ndarray:
    """Raw cross-moment spike coordinates zeta_j = sqrt(n) v_j'b / sqrt(d_j).

    b = Xc'Yc/n with Yc CENTERED but NOT standardized: the coordinate mean
    shift under H1 is v_j'Lambda gamma (absolute units), free of any
    response-scale dilution."""
    n = len(Yc)
    d, V = eig
    k = min(max(int(ktop), 1), len(d))
    b = Xc.T @ Yc / n
    return np.sqrt(n) * (V[:, :k].T @ b) / np.sqrt(d[:k])


def bench_alarm_stat(Xc: np.ndarray, Yc: np.ndarray, eig, bench: Bench,
                     coord_scales: np.ndarray | None = None):
    """Phase-3 gate alarm (deviation D-B0, pre-pass-1 freeze).

    T = max_{j < ktop} |zeta_j| / s_j with s_j the PER-COORDINATE EMPIRICAL
    NULL SCALES s_j = sqrt(mean_null(zeta_j^2)) estimated from the matched
    twins (pass 1) and frozen before any positive arm runs. Rationale: real
    benchmark geometries violate the MP-white-bulk premise of the F12 law
    (family A's smooth bulk, huge spikes), so analytic coordinate variances
    misestimate by geometry-dependent constants; the twin-estimated scales
    absorb every such factor by construction. Before scales are available
    (inside pass 1 itself) s_j = 1 is used and only the pooled mc95 of T
    matters for thresholding; scale estimation consumes the SAME null pool
    (in-sample for size, out-of-sample for every control/positive arm).

    Returns (T, ktop_used)."""
    zeta = raw_z_coords(Xc, Yc, eig, bench.ktop)
    if coord_scales is None:
        s = np.ones(len(zeta))
    else:
        s = np.maximum(np.asarray(coord_scales, float), 1e-12)
    return float(np.max(np.abs(zeta) / s)), int(len(zeta))


def js_asym_from_stats(eig, Xc: np.ndarray, Yc: np.ndarray,
                       K: int | None = None) -> float:
    """Janzing-Schoelkopf-style asymmetry (APPROXIMATE transcription).

    C_R = C - m m'/var(Y) removes the response-explained rank-one component;
    its top eigenvalues interlace with d_i and the relative drops
    (d_i - lam_i)/d_i measure how concentrated the response-explained
    variation is on dominant design directions. Statistic = max drop over
    i < K (K = ktop + 2). Calibrated by permutation upstream."""
    d, V = eig
    n = len(Yc)
    m = Xc.T @ Yc / n
    sy = float(np.mean(Yc ** 2))
    if sy <= 0:
        return 0.0
    w = V.T @ m
    K = K or 4
    K = int(min(K, len(d) - 1))

    def h(lam):
        return float(np.sum(w ** 2 / (sy * (d - lam))))

    drops = []
    for i in range(K):
        lo, hi = d[i + 1], d[i]
        if hi - lo < 1e-14 * max(hi, 1e-300):
            continue
        # h is strictly increasing on (lo, hi) with h(hi^-) -> +inf iff
        # w_i != 0; a root of h(lam) = 1 exists iff h just below hi > 1.
        eps_hi = 1e-9 * (hi - lo)
        if h(hi - eps_hi) <= 1.0:
            continue
        a_, b_ = lo + eps_hi, hi - eps_hi
        for _ in range(80):
            mid = 0.5 * (a_ + b_)
            if h(mid) > 1.0:
                b_ = mid
            else:
                a_ = mid
            if b_ - a_ < 1e-11 * hi:
                break
        lam = 0.5 * (a_ + b_)
        drops.append((d[i] - lam) / d[i])
    return float(max(drops)) if drops else 0.0


def tau_estimators(Xc, Dc, Yc, eig, bench: Bench) -> dict:
    from de_formulas import onatski_select as _on

    d, V = eig

    def tau_joint(A):
        coef, *_ = np.linalg.lstsq(A, Yc, rcond=None)
        return float(coef[0])

    out = {}
    out["tau_ols"] = tau_joint(np.column_stack([Dc, Xc]))
    k = int(max(_on(d), bench.r_inj))
    S = Xc @ V[:, :k]
    out["tau_trim_onatski"] = tau_joint(np.column_stack([Dc, S]))
    lam = 1.0
    n = len(Dc)
    Sc = Xc.T @ Xc / n + lam * np.eye(Xc.shape[1])
    rhs = Xc.T @ Yc / n
    c_vec = Xc.T @ Dc / n
    m = float(Dc @ Dc) / n
    Scinv_c = np.linalg.solve(Sc, c_vec)
    Scinv_rhs = np.linalg.solve(Sc, rhs)
    out["tau_ridge1"] = (
        float(Dc @ Yc) / n - c_vec @ Scinv_rhs
    ) / (m - c_vec @ Scinv_c)
    out["k_trim"] = k
    return out


# ---------------------------------------------------------------------------
# frontier prediction g* (F12 law; mirrors phase2_analysis construction)
# ---------------------------------------------------------------------------


def predicted_g_star(bench: Bench, mc95: float,
                     coord_scales: np.ndarray | None = None,
                     r: int | None = None, seed: int = 0) -> float:
    """F12-law frontier for the benchmark alarm (deviation D-B0 form).

    Per-unit-g mean shift of the raw coordinate zeta_j under H1:
        m_j(g) = g * dir_j * omega_j * sqrt(n) * sqrt(se2 l_j) / sqrt(d_pred)
    with d_pred = bbp_location(l_j, c, se2) and omega the clipped
    min-norm-capture weight at c > 1. Coordinate noise scales are the
    EMPIRICAL twin scales s_j (D-B0), so the standardized shift is
    m_j(g)/s_j and the max-statistic power curve follows the same MC
    construction as phase2_analysis.predicted_frontier_g.
    """
    c, n = bench.c, bench.n
    r = int(r if r is not None else bench.r_inj)
    ktop_eff = max(bench.ktop, r)
    sup = list(range(min(r, ktop_eff)))
    if not sup:
        return float("inf")
    dirv = np.ones(len(sup)) / np.sqrt(len(sup))
    if coord_scales is None:
        s = np.ones(len(sup))
    else:
        s = np.maximum(np.asarray(coord_scales, float)[:len(sup)], 1e-12)
    slope = np.zeros(len(sup))
    for i, j in enumerate(sup):
        lj = float(bench.l_hat[j])
        dj = bbp_location(lj, c, bench.se2)
        omega = float(np.clip(minnorm_capture(np.array([lj]), c)[0], 0, 1)) \
            if c > 1 else 1.0
        slope[i] = (np.sqrt(n) * omega * np.sqrt(bench.se2 * lj) * dirv[i] /
                    (np.sqrt(dj) * s[i]))
    if np.max(slope) <= 0:
        return float("inf")
    rr = np.random.default_rng(seed)
    z0 = np.abs(rr.normal(size=(20000, len(sup))))
    thr_sim = float(np.quantile(z0.max(axis=1), 0.95))
    scale = mc95 / thr_sim
    for g in np.linspace(0.01, 20.0, 400):
        z1 = np.abs(rr.normal(size=(20000, len(sup))) +
                    slope[None, :] * g).max(axis=1)
        if float((z1 * scale > mc95).mean()) >= 0.8:
            return round(float(g), 3)
    return float("inf")


# ---------------------------------------------------------------------------
# arms / cells runner (checkpoint-resume like runners.run_cell)
# ---------------------------------------------------------------------------

ARMS_ORDER = ["null", "perm_null", "pos_half", "pos_1", "pos_2",
              "splithalf", "align_top", "align_weak", "rinj_minus",
              "rinj_plus", "hetero_eps", "m2_null", "m2_pos"]

PASS1_ARMS = ["null", "perm_null"]


def arm_specs():
    import yaml

    cfg = yaml.safe_load((ROOT / "configs" / "benchmarks_frozen.yaml")
                         .read_text())
    return cfg


def skip_arm(bench: Bench, arm: str, spec: dict) -> bool:
    if spec.get("only") and bench.name not in spec["only"]:
        return True
    thr = spec.get("skip_if_r_inj_leq")
    if thr is not None and bench.r_inj <= int(thr):
        return True
    return False


def _atomic_write_parquet(df: pd.DataFrame, path: Path):
    tmp = path.with_suffix(".tmp.parquet")
    df.to_parquet(tmp, index=False)
    tmp.replace(path)


def run_cell(job: dict):
    arm = job["arm"]
    spec = job["spec"]
    raw_path = RAW_DIR / job["config_name"] / f"{arm}.parquet"
    means_path = RAW_DIR / job["config_name"] / f"{arm}_means.npz"
    start_rep = 0
    prev_frames = []
    if raw_path.exists():
        try:
            have = pd.read_parquet(raw_path, columns=["rep"])
            n_have = int(have["rep"].nunique())
            if n_have >= spec["reps"]:
                print(f"[skip] {job['config_name']}/{arm} done ({n_have} reps)",
                      flush=True)
                return arm, 0.0
        except Exception:
            pass
    bench = Bench(job["config_name"])
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    if raw_path.exists():
        try:
            have = pd.read_parquet(raw_path, columns=["rep"])
            n_have = int(have["rep"].nunique())
            if n_have >= spec["reps"]:
                return arm, 0.0
            start_rep = n_have
            prev_frames = [pd.read_parquet(raw_path)]
        except Exception:
            start_rep = 0
    t0 = time.perf_counter()
    frames = list(prev_frames)
    acc_sum = np.zeros(bench.p)
    acc_n = 0
    g_star = job.get("g_star")
    cs = job.get("coord_scales")
    since_flush = 0
    if start_rep and means_path.exists():
        try:
            mz = np.load(means_path)
            acc_sum[: len(mz["mean_bias"])] += (
                np.asarray(mz["mean_bias"], float) * float(mz["n_reps"]))
            acc_n = int(mz["n_reps"])
        except Exception:
            acc_n = 0
    for rep in range(start_rep, spec["reps"]):
        rows, mean_diff = run_rep(bench, bench.name, spec, arm, rep, g_star,
                                  coord_scales=cs)
        frames.append(pd.DataFrame(rows))
        if mean_diff is not None:
            acc_sum += mean_diff
            acc_n += 1
        since_flush += 1
        if since_flush >= 20:
            _atomic_write_parquet(pd.concat(frames, ignore_index=True),
                                  raw_path)
            if acc_n:
                np.savez_compressed(means_path,
                                    mean_bias=(acc_sum / acc_n),
                                    n_reps=acc_n)
            since_flush = 0
    df = pd.concat(frames, ignore_index=True)
    _atomic_write_parquet(df, raw_path)
    if acc_n:
        np.savez_compressed(means_path,
                            mean_bias=(acc_sum / acc_n).astype(np.float64),
                            n_reps=acc_n)
    dt = time.perf_counter() - t0
    print(f"[done] {bench.name}/{arm}: {len(df)} rows (+{spec['reps'] - start_rep}) "
          f"in {dt:.1f}s", flush=True)
    return f"{bench.name}/{arm}", dt


def build_jobs(pass_name: str, freeze: dict | None, workers: int = 3):
    cfg = arm_specs()
    audit = json.loads((BENCH_DIR / "spectral_audit.json").read_text())
    profiles = audit.get("spectral_profiles", audit)
    jobs = []
    for name, ccfg in cfg["configs"].items():
        r_inj = int(profiles[name]["r_inj"])
        for arm in ARMS_ORDER:
            spec = cfg["arms"][arm]
            thr = spec.get("skip_if_r_inj_leq")
            if spec.get("only") and name not in spec["only"]:
                continue
            if thr is not None and r_inj <= int(thr):
                continue
            want = PASS1_ARMS if pass_name == "pass1" else \
                [a for a in ARMS_ORDER if a not in PASS1_ARMS]
            if arm not in want:
                continue
            gs = None
            cs = None
            if pass_name == "pass2" and freeze:
                entry = freeze["configs"][name]
                gs = entry["g_star"]
                cs = np.asarray(entry["coord_scales"], float)
            jobs.append({"config_name": name, "arm": arm, "spec": spec,
                         "g_star": gs, "coord_scales": cs})
    return jobs


def run_jobs(jobs: list[dict], workers: int = 3):
    from multiprocessing import Pool

    if workers <= 1:
        results = [run_cell(j) for j in jobs]
    else:
        with Pool(workers) as pool:
            results = pool.map(run_cell, jobs)
    return results
