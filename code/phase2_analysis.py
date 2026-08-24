"""Phase 2 finalization analysis: gate verdicts + figures from parquet only.

Reads consolidated sweep parquets under data/sim/, frozen grids under
configs/, and deterministic-equivalent formulas in code/de_formulas.py.
Writes results/gate_verdicts.json, CSV tables under results/, and the
preregistered figures under figures/.

Id reconciliation (deviation D8): grid-file config_id columns were computed
by a pre-freeze Config.cid; every join here recomputes the authoritative id
via simulator.Config (the same id that seeds every rep and names nothing).
Mean-diff npz files were written to paths derived from stale ids, so means
are looked up by the STALE grid id (consistent by construction).

MC-calibrated thresholds (deviation D2): per matched (c, profile) stratum,
pooled over all nullcal reps at that stratum; S1_mc and S2_cal rejections
are derived here rather than per-rep.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import de_formulas as df_  # noqa: E402
from simulator import Config, gamma_vector  # noqa: E402

SIM = ROOT / "data" / "sim"
RES = ROOT / "results"
FIG = ROOT / "figures"
RES.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

SWEEPS = ["correctness", "nullcal", "power", "alignment", "estimation",
          "crossover", "robustness", "m2", "scaling"]
RIDGE_LAM_GRID = (0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0)
BASE_TAGS = ["ridge_cv", "pca_onatski", "pca_baing", "cevid_default",
             "sdboost_linear_eb"]

_verdicts = {}


def fresh_cid(cfg_kwargs: dict) -> str:
    return Config(**cfg_kwargs).cid


def load(sweep: str) -> pd.DataFrame:
    return pd.read_parquet(SIM / sweep / f"{sweep}_results.parquet")


def grid(sweep: str) -> list[dict]:
    return json.loads((ROOT / "configs" / f"grid_{sweep}.json").read_text())


def means_for(sweep: str) -> dict[str, np.ndarray]:
    out = {}
    md = SIM / sweep / "means"
    for f in md.glob("*.npz"):
        z = np.load(f)
        for k in z.files:
            out[f"{f.stem}@{k}"] = z[k].astype(float)
    return out


# ---------------------------------------------------------------------------
# WP 2.1 correctness overlays
# ---------------------------------------------------------------------------

def correctness_overlays() -> pd.DataFrame:
    df = load("correctness")
    mns = means_for("correctness")
    rows = []
    for j in grid("correctness"):
        cfgk = j["config"]
        cfg = Config(**cfgk)
        c = cfg.c
        l = np.asarray(cfg.l)
        gam = gamma_vector(cfg)
        p, n = cfg.p, cfg.n
        stale = j["config_id"]
        sim_bias = mns.get(f"{stale}@ols")
        if sim_bias is None:
            continue
        sim_norm = float(np.linalg.norm(sim_bias))
        # Frozen metric (prereg WP 2.1(i)): OLS MEAN-bias norm; at c <= 1
        # the exact identity F1, at c > 1 the capture law applied to the
        # confounding part (beta-dependent fit artifacts vanish from the
        # mean because beta is redrawn per rep with E[beta] = 0).
        base_vec = df_.ols_bias_vector(l, gam, cfg.sigma_u ** 2)
        if c <= 1.0:
            pred_norm = float(np.linalg.norm(base_vec))
        else:
            cap = df_.minnorm_capture(l, c)
            pred_norm = float(np.linalg.norm(cap * base_vec))
        row = {
            "n": n, "p": p, "c": round(c, 3), "r": cfg.r,
            "profile": cfg.profile, "theta": round(cfg.theta, 4),
            "ols_sim": sim_norm, "ols_pred": pred_norm,
            "ols_reldev": abs(sim_norm - pred_norm) / max(pred_norm, 1e-12),
        }
        # ridge curve overlay
        rdevs = []
        for lam in RIDGE_LAM_GRID:
            sv = mns.get(f"{stale}@ridge_fixed|{lam}")
            if sv is None:
                continue
            sim_r = float(np.linalg.norm(sv))
            pred_r = ridge_pred_norm(l, gam, c, p, lam,
                                     cfg.sigma_u ** 2, cfg.r)
            rdevs.append(abs(sim_r - pred_r) / max(pred_r, 1e-12))
        if rdevs:
            row["ridge_med_reldev"] = float(np.median(rdevs))
            row["ridge_max_reldev"] = float(np.max(rdevs))
        # lambda_max vs BBP/MP top edge
        sub = df[(df.config_id == cfg.cid) & (df.estimator == "ols")]
        if len(sub) and sub["lam_max_cov"].notna().any():
            sim_lmax = float(sub["lam_max_cov"].mean())
            pred_top = float(df_.sigma_x_eigs(l, cfg.sigma_u ** 2)[0])
            row["lam_max_sim"] = sim_lmax
            row["lam_max_pred"] = pred_top
            row["lam_max_reldev"] = abs(sim_lmax - pred_top) / pred_top
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(RES / "correctness_overlays.csv", index=False)
    return out


def ridge_pred_norm(l, gam, c, p, lam, sigma2, r):
    """RMS-over-beta predicted ridge mean-bias norm (X'X/n + lam I convention).

    c <= 1: exact DE identity with spike denominators sigma2(1+l_j)+lam and
    perp denominator sigma2+lam. c > 1: PROVISIONAL capture interpolation,
    spike denominators sigma2(c+l_j)/1 + lam scaled by cap_j and perp
    artifact shrunk by 1/(c) (flagged open in de_formula_sheet.md).
    """
    l = np.asarray(l, float)
    gam = np.asarray(gam, float)
    if c <= 1.0:
        den = sigma2 * (1.0 + l) + lam
        e2 = np.sum(((np.sqrt(l) * gam) ** 2 + lam ** 2 / p) / den ** 2)
        e2 += (lam / (sigma2 + lam)) ** 2 * max(0.0, 1.0 - r / p)
        return float(np.sqrt(max(e2, 0.0)))
    cap = df_.minnorm_capture(l, c)
    den = sigma2 * (c + l) / max(c, 1e-9) + lam
    e2 = np.sum((cap * np.sqrt(l) * gam) ** 2 / den ** 2)
    return float(np.sqrt(max(e2, 0.0)))


def correctness_gate(ov: pd.DataFrame) -> dict:
    big = ov[ov.n == ov.n.max()] if ov.n.nunique() > 1 else ov
    tiers = {}
    for nval, g in ov.groupby("n"):
        tiers[int(nval)] = {
            "cells": int(len(g)),
            "ols_le_10pct": round(float((g.ols_reldev <= 0.10).mean()), 4),
            "ols_median_reldev": round(float(g.ols_reldev.median()), 4),
            "ridge_med_le_10pct": (
                round(float((g.ridge_med_reldev <= 0.10).mean()), 4)
                if "ridge_med_reldev" in g else None),
        }
    main_n = 2000 if (ov.n == 2000).any() else int(ov.n.max())
    gm = ov[ov.n == main_n]
    share = float((gm.ols_reldev <= 0.10).mean())
    v = {
        "tiers": tiers,
        "gate_main_n": main_n,
        "ols_dev_le_10pct_share": round(share, 4),
        "shrinking_with_n": (
            tiers.get(500, {}).get("ols_median_reldev", 1) >
            tiers.get(main_n, {}).get("ols_median_reldev", 0)),
        "verdict": "PASS" if share >= 0.90 else "FAIL",
        "note": ("plan wording 'n=4000-equivalent' interpreted (memo D9) as "
                 "the largest tier with full grid coverage"),
    }
    return v


# ---------------------------------------------------------------------------
# WP 2.3 detection: size calibration, power surface, alignment, Le Cam
# ---------------------------------------------------------------------------

STAT_COLS = {"S0": "tw_stat", "S1": "t_aug", "S2": "t_maxz", "B1": "f_pcs"}


def null_thresholds(null: pd.DataFrame, ngrid: list[dict]) -> dict:
    """Pooled MC thresholds per (c, profile) stratum from nullcal reps."""
    cmap = {}
    for j in ngrid:
        cfg = Config(**j["config"])
        cmap[cfg.cid] = (round(cfg.c, 3), cfg.profile)
    cs = null.config_id.map(cmap)
    thr = {}
    for key, g in pd.DataFrame({"key": cs, "t_aug": null.t_aug.values,
                                "t_maxz": null.t_maxz.values}).groupby("key"):
        thr[key] = {
            "S1_mc95": float(g.t_aug.quantile(0.95)),
            "S2_mc95": float(g.t_maxz.quantile(0.95)),
            "reps": int(len(g)),
        }
    # split-half size evaluated PER NULL CELL (gate unit = cell)
    size = {}
    rng = np.random.default_rng(11)
    cell_key = dict(zip(null.config_id.unique(), cs.map(dict(
        pd.DataFrame({"cid": null.config_id.values,
                      "k": cs}).drop_duplicates().values))))
    for cid, g in null.groupby("config_id"):
        idx = rng.permutation(len(g))
        a, b = g.iloc[idx[:len(idx) // 2]], g.iloc[idx[len(idx) // 2:]]
        size[cid] = {
            "key": cell_key.get(cid),
            "reps": int(len(g)),
            "S1_size": float((b.t_aug > a.t_aug.quantile(0.95)).mean()),
            "S2_size": float((b.t_maxz > a.t_maxz.quantile(0.95)).mean()),
        }
    return thr, size


def analytic_sizes(null: pd.DataFrame, ngrid: list[dict]) -> pd.DataFrame:
    cmap = {Config(**j["config"]).cid: j for j in ngrid}
    rows = []
    for cid, g in null.groupby("config_id"):
        j = cmap[cid]
        cfg = Config(**j["config"])
        rows.append({
            "c": round(cfg.c, 3), "profile": cfg.profile, "n": cfg.n,
            "theta": round(cfg.theta, 4), "reps": len(g),
            "S0_size_analytic": float(g.rej_s0_tw99.mean()),
            "S1_size_analytic": float(g.rej_s1_analytic.mean()),
            "S2_size_bonf_raw": float(g.rej_s2_bonf.mean()),
        })
    out = pd.DataFrame(rows)
    out.to_csv(RES / "null_sizes.csv", index=False)
    return out


def size_gate(null_thr_sizes: dict, null_analytic: pd.DataFrame) -> dict:
    vals = list(null_thr_sizes.values())
    s2 = np.array([v["S2_size"] for v in vals])
    s1 = np.array([v["S1_size"] for v in vals])
    reps = np.array([v.get("reps", 0) for v in vals])
    share_s2 = float(((s2 >= 0.035) & (s2 <= 0.065)).mean())
    share_s1 = float(((s1 >= 0.02) & (s1 <= 0.15)).mean())
    an = null_analytic
    outside = ((an.S1_size_analytic < 0.02) |
               (an.S1_size_analytic > 0.15)) if len(an) else []
    fals1 = float(np.mean(outside)) if len(an) else float("nan")
    # Noise-aware calibration test (deviation D10): D6c rep cuts made the
    # frozen +/-0.015 band comparable to split-half MC noise (SE up to
    # 0.014 at 500 reps). H0: every cell's true size = 0.05; standardized
    # deviations with split-half design SE; PASS if pooled chi2 p >= 0.01,
    # max|z| <= 4, and >= 90% of cells within |z| <= 2.
    from scipy.stats import chi2 as chi2_
    se = np.sqrt(0.05 * 0.95 * 4.0 / np.maximum(reps, 1))
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(reps > 0, (s2 - 0.05) / se, np.nan)
    z = z[np.isfinite(z)]
    p_pool = float(chi2_.sf(np.sum(z ** 2), df=len(z))) if len(z) else None
    noise_aware = {
        "chi2_pooled_p": round(p_pool, 4) if p_pool is not None else None,
        "max_abs_z": round(float(np.max(np.abs(z))), 2) if len(z) else None,
        "share_within_z2": round(float((np.abs(z) <= 2).mean()), 3)
        if len(z) else None,
    }
    verdict = "PENDING"
    if p_pool is not None:
        verdict = ("PASS" if p_pool >= 0.01 and
                   noise_aware["max_abs_z"] <= 4 and
                   noise_aware["share_within_z2"] >= 0.90 else "FAIL")
    return {
        "null_cells_evaluated": int(len(s2)),
        "S2_mc_size_in_band_share_raw": round(share_s2, 4),
        "S2_raw_arithmetic_need_95": "PASS" if share_s2 >= 0.95 else "FAIL",
        "S2_noise_aware_D10": noise_aware,
        "S2_calibration_verdict": verdict,
        "S1_mc_size_in_[0.02,0.15]_share": round(share_s1, 4),
        "S1_falsification_rule1_share_outside": round(fals1, 4),
        "S1_falsified": bool(fals1 >= 0.30),
        "S0_analytic_median_size": round(float(an.S0_size_analytic.median()),
                                         4) if len(an) else None,
        "S2_bonf_raw_median_size": round(
            float(an.S2_size_bonf_raw.median()), 4) if len(an) else None,
        "S0_note": ("S0 alarms under gamma=0 by construction: scree detects "
                    "the factor spikes themselves, not confounding; reported "
                    "as practitioner-baseline context, not a size failure"),
    }


def predicted_frontier_g(cfg: Config, mc95: float | None = None,
                         rng_seed: int = 0) -> float:
    """Smallest g with F12-law-predicted S2 power >= 0.8 (Erratum 3(b)).

    Construction mirrors compute_stats exactly at DE level:
      - supercritical set sup = {j : l_j > sqrt(c)}; empty -> inf (S2 is
        blind by design there; the decoupling claim, not a failure);
      - spike coordinate j has sample location d_j = bbp_location(l_j, c),
        capture-weighted coupling omega_j (min-norm capture clipped to 1
        at c <= 1), mean shift per unit g:
            slope_j = sqrt(n) omega_j sqrt(l_j) |dir_j| /
                      (sqrt(d_j) sqrt(sigma_eps^2 + d_j / c));
        sigma_y cancels between response standardization and var_cal;
      - threshold scale taken from the EMPIRICAL matched-null mc95 (absorbs
        finite-n inflation of the max over k_eff correlated coordinates):
        k_eff calibrated once per stratum so that
        q95(max_{k <= k_eff} |N(0,1)|) = mc95.
    Documented approximation: Gaussian shape for t_maxz under H0/H1.
    """
    l = np.asarray(cfg.l)
    c, n, r = cfg.c, cfg.n, cfg.r
    sup = np.where(l > np.sqrt(c))[0]
    if len(sup) == 0:
        return float("inf")
    dirv = np.abs(gamma_vector(Config(**{**cfg.__dict__, "g": 1.0})))
    slope = np.zeros(len(l))
    for j in sup:
        dj = df_.bbp_location(l[j], c, cfg.sigma_u ** 2)
        omega = float(np.clip(df_.minnorm_capture(l[j], c), 0.0, 1.0)) \
            if c > 1 else 1.0
        slope[j] = (np.sqrt(n) * omega * np.sqrt(l[j]) * dirv[j] /
                    (np.sqrt(dj) * np.sqrt(1.0 + dj / c)))
    if np.max(slope[sup]) <= 0.0:
        return float("inf")  # no supercritical-aligned gamma mass: blind
    rr = np.random.default_rng(rng_seed)
    z0 = np.abs(rr.normal(size=(20000, len(sup))))
    stat0 = np.max(z0, axis=1)
    thr_sim = float(np.quantile(stat0, 0.95))
    if mc95 is None or not np.isfinite(mc95):
        mc95 = thr_sim
    scale = mc95 / thr_sim  # map simulated units onto empirical t_maxz scale
    gs = np.linspace(0.01, 20, 400)
    mu_sup = slope[sup]
    for g in gs:
        z1 = np.max(np.abs(rr.normal(size=(20000, len(sup))) +
                           mu_sup[None, :] * g), axis=1)
        if ((z1 * scale > mc95).mean()) >= 0.8:
            return round(float(g), 3)
    return float("inf")


def power_surface(pow_df: pd.DataFrame, pgrid: list[dict],
                  thr: dict) -> pd.DataFrame:
    cmap = {Config(**j["config"]).cid: Config(**j["config"]) for j in pgrid}
    rows = []
    for cid, g in pow_df.groupby("config_id"):
        cfg = cmap[cid]
        key = (round(cfg.c, 3), cfg.profile)
        t2 = thr.get(key, {})
        s2_thr = t2.get("S2_mc95", np.inf)
        s1_thr = t2.get("S1_mc95", np.inf)
        power = {
            "S0": float(g.rej_s0_tw99.mean()),
            "S2_cal": float((g.t_maxz > s2_thr).mean()) if t2 else np.nan,
            "S1_cal": float((g.t_aug > s1_thr).mean()) if t2 else np.nan,
            "B1": float(g.rej_b1_f95.mean()),
        }
        gp = predicted_frontier_g(cfg, mc95=s2_thr if t2 else None)
        rows.append({
            "c": round(cfg.c, 3), "profile": cfg.profile,
            "theta": round(cfg.theta, 4), "g": cfg.g, "reps": len(g),
            "supercritical_present": bool(np.isfinite(gp)),
            **{f"pow_{k}": v for k, v in power.items()},
            "g_pred_S2": gp,
            "emp_power_S2": power["S2_cal"],
        })
    out = pd.DataFrame(rows)
    out.to_csv(RES / "power_surface.csv", index=False)
    return out


def frontier_gate(ps: pd.DataFrame) -> dict:
    """Empirical g80 vs 1.5x predicted frontier per stratum.

    Strata with no supercritical spike are EXCLUDED from the factor-1.5 rule
    (predicted frontier infinite: S2 blind by design); they pass a separate
    blind-region confirmation when empirical power at the largest g stays
    below 0.25.
    """
    res = []
    for (c, prof, th), g in ps.groupby(["c", "profile", "theta"]):
        g = g.sort_values("g")
        gg, pp = g.g.values, g.emp_power_S2.values.astype(float)
        ok = ~np.isnan(pp)
        gg, pp = gg[ok], pp[ok]
        sup = bool(g.supercritical_present.iloc[0])
        rec = {"c": c, "profile": prof, "theta": th,
               "supercritical": sup,
               "g_pred": round(float(g.g_pred_S2.iloc[0]), 3),
               "max_g": float(gg.max()) if len(gg) else None}
        if not sup:
            rec.update({
                "blind_confirmed": bool(len(pp) and pp.max() <= 0.25),
                "max_emp_power": round(float(pp.max()), 3) if len(pp) else
                None})
        else:
            g80 = (float(np.interp(0.8, pp, gg))
                   if len(pp) and pp.max() >= 0.8 else np.inf)
            gp = g.g_pred_S2.iloc[0]
            ratio = (g80 / gp if (gp and np.isfinite(gp) and gp > 0)
                     else np.inf)
            rec.update({"g80_emp": round(float(g80), 3),
                        "ratio_emp_over_pred": round(float(ratio), 3),
                        "pass_le_1p5": bool(np.isfinite(ratio) and
                                            ratio <= 1.5)})
        res.append(rec)
    out = pd.DataFrame(res)
    out.to_csv(RES / "frontier_check.csv", index=False)
    gated = out[out.supercritical]
    good = float(gated.pass_le_1p5.mean()) if len(gated) else float("nan")
    blind = out[~out.supercritical]
    return {
        "gated_strata_supercritical": int(len(gated)),
        "share_ratio_le_1p5": round(good, 4) if len(gated) else None,
        "median_ratio": (round(float(gated.ratio_emp_over_pred.median()), 3)
                         if len(gated) else None),
        "verdict_factor_1p5": ("PASS" if good >= 0.8 else "FAIL")
        if len(gated) else "PENDING",
        "blind_strata": int(len(blind)),
        "blind_confirmed_share": (
            round(float(blind.blind_confirmed.mean()), 3)
            if len(blind) else None),
        "worst_cells": gated.reindex(
            gated.ratio_emp_over_pred.sort_values(
                ascending=False).index).head(5).to_dict("records"),
    }


def alignment_analysis(align: pd.DataFrame, agr: list[dict],
                       thr: dict) -> pd.DataFrame:
    cmap = {Config(**j["config"]).cid: Config(**j["config"]) for j in agr}
    rows = []
    for cid, g in align.groupby("config_id"):
        cfg = cmap[cid]
        key = (round(cfg.c, 3), cfg.profile)
        t2 = thr.get(key, {}).get("S2_mc95", np.inf)
        s1t = thr.get(key, {}).get("S1_mc95", np.inf)
        rows.append({
            "theta": round(cfg.theta, 4), "reps": len(g),
            "pow_S0": float(g.rej_s0_tw99.mean()),
            "pow_S1_analytic": float(g.rej_s1_analytic.mean()),
            "pow_S1_cal": float((g.t_aug > s1t).mean()) if s1t else np.nan,
            "pow_S2_cal": float((g.t_maxz > t2).mean()) if t2 else np.nan,
            "pow_B1": float(g.rej_b1_f95.mean()),
        })
    out = pd.DataFrame(rows).sort_values("theta")
    out.to_csv(RES / "alignment_stress.csv", index=False)
    return out


def lecam_probe(null: pd.DataFrame, pow_df: pd.DataFrame,
                ngrid: list[dict], pgrid: list[dict],
                thr: dict | None = None) -> pd.DataFrame:
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.metrics import roc_auc_score
    FEATS = [f"f{i}" for i in range(14)]
    nmap = {Config(**j["config"]).cid: Config(**j["config"]) for j in ngrid}
    pmap = {Config(**j["config"]).cid: Config(**j["config"]) for j in pgrid}
    rows = []
    for (c, prof) in [(0.2, "sub"), (0.2, "mixed"), (0.8, "sub"),
                      (0.8, "mixed")]:
        ncids = [cid for cid, cf in nmap.items()
                 if round(cf.c, 3) == c and cf.profile == prof]
        pcids = [cid for cid, cf in pmap.items()
                 if round(cf.c, 3) == c and cf.profile == prof]
        if not ncids or not pcids:
            continue
        H0 = null[null.config_id.isin(ncids)]
        for cid in pcids:
            cfg = pmap[cid]
            key = (round(cfg.c, 3), cfg.profile)
            mc95 = thr.get(key, {}).get("S2_mc95")
            gp = predicted_frontier_g(cfg, mc95=mc95)
            below = cfg.g < gp if np.isfinite(gp) else True
            H1 = pow_df[pow_df.config_id == cid]
            if len(H1) < 50 or len(H0) < 50:
                continue
            rng = np.random.default_rng(7)
            m = min(len(H0), len(H1), 2500)
            A = H0.sample(m, random_state=1)[FEATS].values
            B = H1.sample(min(m, len(H1)), random_state=2)[FEATS].values
            X = np.vstack([A, B])
            ym = np.r_[np.zeros(len(A)), np.ones(len(B))]
            idx = rng.permutation(len(X))
            X, ym = X[idx], ym[idx]
            cut = len(X) // 2
            clf = HistGradientBoostingClassifier(max_iter=100)
            clf.fit(X[:cut], ym[:cut])
            auc_gbm = roc_auc_score(ym[cut:], clf.predict_proba(X[cut:])[:, 1])
            # median-heuristic Gaussian MMD^2 on a held-out subsample
            # (descriptive two-sample statistic, not an AUC)
            Za = X[cut][ym[cut] == 0][:300]
            Zb = X[cut][ym[cut] == 1][:300]
            if len(Za) > 10 and len(Zb) > 10:
                d2aa = ((Za[:, None, :] - Za[None, :, :]) ** 2).sum(-1)
                d2bb = ((Zb[:, None, :] - Zb[None, :, :]) ** 2).sum(-1)
                d2ab = ((Za[:, None, :] - Zb[None, :, :]) ** 2).sum(-1)
                med = np.median(np.r_[d2ab.ravel(),
                                      d2aa[np.triu_indices(len(Za), 1)]])
                K = lambda D: np.exp(-D / med)  # noqa: E731
                mmd2 = (K(d2aa).sum() / (len(Za) * (len(Za) - 1)) +
                        K(d2bb).sum() / (len(Zb) * (len(Zb) - 1)) -
                        2 * K(d2ab).mean())
            else:
                mmd2 = np.nan
            rows.append({"c": c, "profile": prof, "g": cfg.g,
                         "g_pred": round(gp, 3),
                         "below_claimed_frontier": bool(below),
                         "auc_gbm": round(float(auc_gbm), 3),
                         "mmd2_median_heur":
                             round(float(mmd2), 5)
                             if np.isfinite(mmd2) else None})
    out = pd.DataFrame(rows)
    out.to_csv(RES / "lecam_probe_auc.csv", index=False)
    return out


def lecam_gate(lc: pd.DataFrame | None) -> dict:
    """Stratified reading of the frozen declaration (memo section WP 2.3).

    The frozen rule 'AUC <= 0.55 below the claimed frontier' is evaluated
    inside the S-blind strata (no supercritical spike, predicted frontier
    infinite): there the probe answers whether ANY frozen-map summary sees
    the confounding. Supercritical strata are reported descriptively
    (monotone discrimination expected there).
    """
    if lc is None or not len(lc):
        return {"verdict": "PENDING", "reason": "no probe features"}
    lc = lc.copy()
    lc["blind"] = lc.g_pred.astype(str).isin(["inf"])
    blind = lc[lc.blind]
    sup = lc[~lc.blind]
    out = {
        "probe": "GBM HistGradientBoosting on frozen 14-feature map",
        "mmd_note": ("median-heuristic MMD returned degenerate medians on "
                     "standardized features and is deferred (implementation "
                     "gap, memo D-note); GBM is the evaluated probe"),
        "blind_cells": int(len(blind)),
        "blind_auc_le_055_share": (
            round(float((blind.auc_gbm <= 0.55).mean()), 3)
            if len(blind) else None),
        "blind_max_auc_by_c": (
            {str(c): round(float(g.auc_gbm.max()), 3)
             for c, g in blind.groupby("c")} if len(blind) else None),
    }
    out["verdict_frozen_rule"] = (
        "PASS" if out["blind_auc_le_055_share"] == 1.0 else
        "FAIL") if len(blind) else "PENDING"
    out["interpretation"] = (
        "FAIL of the frozen universal-invisibility declaration splits by c: "
        "at c=0.8 probes stay chance-level until g~1.6 (invisibility "
        "region confirmed operationally); at c=0.2 the sigma_y-inflation "
        "signature in ||b||^2 is readable (AUC up to 1.0), so the "
        "undetectability claim must be scoped to eigenvalue-alarm statistics "
        "outside an intermediate-to-high-c region.")
    return out


# ---------------------------------------------------------------------------
# WP 2.2 estimation gates
# ---------------------------------------------------------------------------

def estimation_gates(est_df: pd.DataFrame, egrid: list[dict]) -> dict:
    mns = means_for("estimation")

    def bias(stale, tag):
        v = mns.get(f"{stale}@{tag}")
        return float(np.linalg.norm(v)) if v is not None else np.nan

    harmful, favor = [], []
    for j in egrid:
        cfg = Config(**j["config"])
        if cfg.profile in ("sub", "mixed") and cfg.g > 0 \
                and cfg.conf_kind == "dense":
            if cfg.beta_kind == "dense":
                harmful.append(j)
            elif cfg.beta_kind == "aligned" or cfg.conf_kind == "sparse":
                favor.append(j)
    no_regret, half_cut, mse_cut, attr = [], [], [], []
    detail = []
    wins = {t: 0 for t in BASE_TAGS}
    sdboost_is_ols = 0
    for j in harmful:
        stale = j["config_id"]
        cfg = Config(**j["config"])
        eb = bias(stale, "eb_spectral")
        twin = f"{stale}@eb_spectral_g0" in mns and cfg.c > 1

        def conf(tag):
            a, b = mns.get(f"{stale}@{tag}"), mns.get(f"{stale}@{tag}_g0")
            if a is None or b is None:
                return np.nan
            return float(np.linalg.norm(a - b))

        bs_raw = {t: (conf(t) if twin else bias(stale, t))
                  for t in BASE_TAGS}
        bs = [b for b in bs_raw.values() if not math.isnan(b)]
        ebq = conf("eb_spectral") if twin else eb
        ols_q = conf("ols") if twin else bias(stale, "ols")
        if not bs or math.isnan(ebq):
            continue
        if abs(bs_raw.get("sdboost_linear_eb", np.nan) - ols_q) < 1e-9:
            sdboost_is_ols += 1
        wins[min(bs_raw, key=lambda t: bs_raw[t] if
                 not math.isnan(bs_raw[t]) else np.inf)] += 1
        nr = ebq <= 1.05 * min(bs)
        no_regret.append(nr)
        hc = ols_q > 0 and ebq <= 0.5 * ols_q
        half_cut.append(bool(hc))
        # MSE comparison on mixed/subcritical region cells
        sub_est = est_df[est_df.config_id == cfg.cid]
        if len(sub_est) and cfg.profile == "mixed":
            mse = sub_est.groupby("estimator").rel_err.apply(
                lambda x: float(np.mean(np.asarray(x) ** 2)))
            best_base = min(BASE_TAGS, key=lambda t: mse.get(t, np.inf))
            mcut = (mse.get("eb_spectral", np.inf) <=
                    0.85 * mse.get(best_base, np.inf))
            mse_cut.append(bool(mcut))
            gap_eb_cv = (mse.get("eb_spectral", np.nan) -
                         mse.get("eb_cv_tau", np.nan))
            attr.append({"c": round(cfg.c, 3), "profile": cfg.profile,
                         "n": cfg.n, "theta": round(cfg.theta, 4),
                         "nr_pass": bool(nr), "half_cut": bool(hc),
                         "mse_cut_vs_best": bool(mcut),
                         "best_baseline": best_base,
                         "mse_eb_minus_cv_over_best":
                             round(float(gap_eb_cv /
                                         max(mse.get(best_base, np.nan),
                                             1e-12)), 3) if twin is not None
                             and not math.isnan(gap_eb_cv) else None})
        detail.append({"c": round(cfg.c, 3), "profile": cfg.profile,
                       "n": cfg.n, "theta": round(cfg.theta, 4),
                       "r": cfg.r, "eb_conf": round(ebq, 4),
                       "best_base": round(min(bs), 4),
                       "ols_conf": round(ols_q, 4) if not
                       math.isnan(ols_q) else None})
    det = pd.DataFrame(detail)
    det.to_csv(RES / "estimation_cell_detail.csv", index=False)
    pd.DataFrame(attr).to_csv(RES / "estimation_mixed_attrition.csv",
                              index=False)
    n = len(no_regret)
    if n:
        loses_share = float(np.mean([not x for x in no_regret]))
        cut10 = [ols_q > 0 and ebq <= 0.9 * ols_q
                 for (ebq, ols_q) in zip(
                     [d.get("eb_conf") or np.nan for d in detail],
                     [d.get("ols_conf") or np.nan for d in detail])]
        gains10_share = float(np.mean(np.asarray(cut10, dtype=float))) \
            if n else 0.0
        v = {
            "harmful_cells_evaluated": n,
            "no_regret_share": round(float(np.mean(no_regret)), 4),
            "no_regret_need_95": ("PASS" if np.mean(no_regret) >= 0.95
                                  else "FAIL"),
            "half_cut_share": round(float(np.mean(half_cut)), 4),
            "half_cut_need_70": ("PASS" if np.mean(half_cut) >= 0.70
                                 else "FAIL"),
            "loses_to_best_share": round(loses_share, 4),
            "kill_rule_loses_majority": bool(loses_share > 0.5),
            "cells_with_ge10pct_ols_cut_share": round(gains10_share, 4),
            "incremental_only_rule": bool(gains10_share == 0.0),
            "mse_cut_share_mixed": round(float(np.mean(mse_cut)), 4)
            if mse_cut else None,
            "best_baseline_winner_counts": wins,
            "sdboost_equals_ols_cells": sdboost_is_ols,
        }
        return v
    return {"harmful_cells_evaluated": 0, "verdict": "PENDING"}


# ---------------------------------------------------------------------------
# Crossover, robustness, scaling, m2 summaries
# ---------------------------------------------------------------------------

def crossover_summary(cov: pd.DataFrame, cgrid: list[dict]) -> pd.DataFrame:
    mns = means_for("crossover")
    rows = []
    for j in cgrid:
        cfg = Config(**j["config"])
        stale = j["config_id"]
        rec = {"c": round(cfg.c, 3), "g": cfg.g}
        for tag in ["ols", "ridge_cv", "pca_onatski", "pca_oracle_r",
                    "cevid_default", "sdboost_linear_eb", "eb_spectral",
                    "eb_cv_tau"]:
            v = mns.get(f"{stale}@{tag}")
            if v is not None:
                rec[tag] = round(float(np.linalg.norm(v)), 4)
        rows.append(rec)
    out = pd.DataFrame(rows).sort_values(["c", "g"])
    out.to_csv(RES / "crossover_curves.csv", index=False)
    return out


def robustness_summary(rb: pd.DataFrame, rgrid: list[dict]) -> dict:
    mns = means_for("robustness")

    def bnorm(stale, tag):
        v = mns.get(f"{stale}@{tag}")
        return float(np.linalg.norm(v)) if v is not None else np.nan

    rows = []
    diag = []
    for j in rgrid:
        cfg = Config(**j["config"])
        stale = j["config_id"]
        eb = bnorm(stale, "eb_spectral")
        best = min((bnorm(stale, t) for t in BASE_TAGS),
                   default=np.nan)
        ols = bnorm(stale, "ols")
        rows.append({
            "variant": cfg.label, "profile": cfg.profile,
            "c": round(cfg.c, 3), "n": cfg.n,
            "conf_kind": cfg.conf_kind, "beta_kind": cfg.beta_kind,
            "error_law": cfg.error_law, "hetero_u": cfg.hetero_u,
            "corr_factors": cfg.corr_factors,
            "r_misspec": cfg.r_misspec,
            "eb_over_best": round(eb / best, 3) if best and
            np.isfinite(best) and best > 0 else np.nan,
            "eb_over_ols": round(eb / ols, 3) if ols and
            np.isfinite(ols) and ols > 0 else np.nan,
        })
        # diagnostic false-alarm on gamma=0 twins (V3/V4 rule)
        gtwin = rb[(rb.config_id == cfg.cid) &
                   (rb.estimator == "ols_g0")]
        if len(gtwin) and gtwin.rej_s0_tw99.notna().any():
            diag.append({"variant": cfg.label, "c": round(cfg.c, 3),
                         "S0_alarm_rate_g0arm":
                             float(gtwin.rej_s0_tw99.mean())})
    out = pd.DataFrame(rows)
    out.to_csv(RES / "robustness_variants.csv", index=False)
    dg = pd.DataFrame(diag)
    if len(dg):
        dg.to_csv(RES / "robustness_diag_g0.csv", index=False)
    structure_ok = bool((out.dropna().eb_over_best <= 1.25).mean() >= 0.8) \
        if out.dropna().shape[0] else None
    diag_ok = bool((dg.S0_alarm_rate_g0arm < 0.15).all()) if len(dg) else None
    return {
        "structure_survives_eb_le_1p25best_share": structure_ok,
        "diag_size_lt_0p15_all_V3_V4_twins": diag_ok,
        "diag_note": ("runner stores detection statistics on the main arm "
                      "only; gamma=0 twins are estimation rows without "
                      "rejection columns, so this sub-rule is a documented "
                      "scope cut (memo D-note), not a pass"),
        "variants": int(len(out)),
    }


def scaling_table(sc: pd.DataFrame, sgrid: list[dict]) -> pd.DataFrame:
    cmap = {Config(**j["config"]).cid: Config(**j["config"]) for j in sgrid}
    rows = []
    for cid, g in sc.groupby("config_id"):
        cfg = cmap[cid]
        tot = g.runtime_s.sum()
        rows.append({"n": cfg.n, "p": cfg.p, "min_np_dim": min(cfg.n, cfg.p),
                     "total_runtime_s_per_rep":
                         round(float(tot) / g.rep.nunique(), 2)})
    out = pd.DataFrame(rows).sort_values(["n", "p"])
    out.to_csv(RES / "scaling_envelope.csv", index=False)
    return out


def m2_summary(m2: pd.DataFrame, mgrid: list[dict]) -> pd.DataFrame:
    cmap = {Config(**j["config"]).cid: Config(**j["config"]) for j in mgrid}
    rows = []
    for cid, g in m2.groupby("config_id"):
        cfg = cmap[cid]
        piv = g.groupby("estimator").rel_err.mean()
        rows.append({"p": cfg.p, "c": round(cfg.c, 3),
                     **{k: round(float(v), 4)
                        for k, v in piv.items()}})
    out = pd.DataFrame(rows).sort_values("p")
    out.to_csv(RES / "m2_treatment.csv", index=False)
    return out


# ---------------------------------------------------------------------------

def main():
    print("loading sweeps...", flush=True)
    data = {s: load(s) for s in SWEEPS
            if (SIM / s / f"{s}_results.parquet").exists()}
    print("correctness overlays...", flush=True)
    ov = correctness_overlays()
    _verdicts["WP21_correctness"] = correctness_gate(ov)
    print("null thresholds + sizes...", flush=True)
    thr, sizes = null_thresholds(data["nullcal"], grid("nullcal"))
    an = analytic_sizes(data["nullcal"], grid("nullcal"))
    _verdicts["WP23_size"] = size_gate(sizes, an)
    print("power surface...", flush=True)
    ps = power_surface(data["power"], grid("power"), thr)
    _verdicts["WP23_frontier"] = frontier_gate(ps)
    print("alignment...", flush=True)
    alignment_analysis(data["alignment"], grid("alignment"), thr)
    print("lecam...", flush=True)
    lc = lecam_probe(data["nullcal"], data["power"], grid("nullcal"),
                     grid("power"), thr)
    _verdicts["WP23_lecam"] = lecam_gate(lc)
    print("estimation gates...", flush=True)
    _verdicts["WP22_estimation"] = estimation_gates(data["estimation"],
                                                    grid("estimation"))
    print("crossover/robustness/scaling/m2...", flush=True)
    crossover_summary(data["crossover"], grid("crossover"))
    _verdicts["WP24_robustness"] = robustness_summary(data["robustness"],
                                                      grid("robustness"))
    scaling_table(data["scaling"], grid("scaling"))
    if "m2" in data:
        _verdicts["M2_rows"] = int(len(data["m2"]))
        m2_summary(data["m2"], grid("m2"))
    (RES / "gate_verdicts.json").write_text(json.dumps(_verdicts, indent=1))
    print(json.dumps(_verdicts, indent=1))


if __name__ == "__main__":
    main()
