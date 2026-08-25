"""SCF Phase 3 analysis (WP 3.3): mechanical gate evaluation from the frozen
thresholds. Consumes data/benchmarks/raw/*/*.parquet +
results/benchmark_freeze.json; produces

  results/benchmark_arms.csv        per config x arm summary
  results/calibration_informality.csv  head-to-head method table
  results/bench_pf_tables.csv       PF-1/PF-3 verdict inputs
  results/bench_g4_verdict.json     mechanical G4 conditions
  figures/benchmark_frontier_check.pdf

No thresholds are estimated here: everything frozen in benchmark_freeze.json
is applied as-is (deviation D-B0 layer).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import benchmarks as B  # noqa: E402

RES = ROOT / "results"
FIG = ROOT / "figures"


def load_all() -> pd.DataFrame:
    frames = []
    for f in sorted(B.RAW_DIR.glob("*/*.parquet")):
        df = pd.read_parquet(f)
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    # Uniform standardized statistic: pass-1 rows carry RAW-scale T_bench
    # (coord scales did not exist yet); pass-2 rows carry standardized ones.
    # Recompute from the zeta columns everywhere so every arm is on the
    # same footing (deviation D-B1 follow-up).
    def std_T(row):
        cn = row["config_name"]
        s = np.asarray(freeze["configs"][cn]["coord_scales"], float)
        k = min(len(s), sum(c.startswith("zeta") for c in df.columns))
        z = np.array([row.get(f"zeta{j}", np.nan) for j in range(k)])
        if np.isnan(z).any():
            return row["T_bench"]
        return float(np.max(np.abs(z) / np.maximum(s[:k], 1e-12)))
    global freeze
    freeze = json.loads((RES / "benchmark_freeze.json").read_text())
    df["T_bench"] = df.apply(std_T, axis=1)
    return df


def main():
    global freeze
    freeze = json.loads((RES / "benchmark_freeze.json").read_text())
    cfg_spec = B.arm_specs()
    df = load_all()
    df["rej_bench"] = [
        t > freeze["configs"][cn]["mc95_S2"]
        for t, cn in zip(df["T_bench"], df["config_name"])]
    df["rej_ucm"] = [
        r > freeze["configs"][cn]["ucm_q95_perm"]
        for r, cn in zip(df["ucm_rho"], df["config_name"])]
    df["rej_js"] = [
        r > freeze["configs"][cn]["js_q95_perm"]
        for r, cn in zip(df["js_asym"], df["config_name"])]

    # ---------------- per-arm summary ---------------------------------
    rows = []
    for (cn, arm), g in df.groupby(["config_name", "arm"]):
        ent = freeze["configs"][cn]
        rec = {
            "config_name": cn, "arm": arm, "n_reps": int(len(g)),
            "size_power_bench": float(g["rej_bench"].mean()),
            "size_power_s0": float(g["rej_s0_tw99"].mean()),
            "size_power_b1": float(g["rej_b1_f95"].mean()),
            "size_power_ucm": float(g["rej_ucm"].mean()),
            "size_power_js": float(g["rej_js"].mean()),
            "T_mean": float(g["T_bench"].mean()),
            "ucm_rho_mean": float(g["ucm_rho"].mean()),
            "js_asym_mean": float(g["js_asym"].mean()),
        }
        for c in ("tau_ols", "tau_trim_onatski", "tau_ridge1"):
            if c in g.columns and g[c].notna().any():
                rec[f"{c}_mean"] = float(g[c].mean())
                rec[f"{c}_err_mean"] = float((g[c] - 1.0).abs().mean())
                rec[f"{c}_sign_ok"] = float((np.sign(g[c]) == 1).mean())
        if "rel_bias_dir" in g.columns and g["rel_bias_dir"].notna().any():
            rec["rel_err_ols_mean"] = float(g["rel_bias_dir"].mean())
        rows.append(rec)
    arms = pd.DataFrame(rows)
    RES.mkdir(exist_ok=True)
    arms.to_csv(RES / "benchmark_arms.csv", index=False)

    # ---------------- PF-1: frontier straddle --------------------------
    pf_rows = []
    for cn, ent in freeze["configs"].items():
        gs = ent["g_star"]

        def pw(arm):
            sel = arms[(arms.config_name == cn) & (arms.arm == arm)]
            return float(sel["size_power_bench"].iloc[0]) if len(sel) \
                else np.nan

        def sz(arm):
            return pw(arm)
        rec = {"config_name": cn, "g_star": gs,
               "size_null_insamp": sz("null"),
               "size_perm": sz("perm_null"),
               "size_splithalf": sz("splithalf"),
               "pow_half": pw("pos_half"), "pow_1": pw("pos_1"),
               "pow_2": pw("pos_2")}
        ctrl_sizes = [v for v in (rec["size_perm"],
                                  rec["size_splithalf"])
                      if np.isfinite(v)]
        rec["max_ctrl_size"] = max(ctrl_sizes) if ctrl_sizes else np.nan
        ok_size = (len(ctrl_sizes) > 0 and 0.02 <= rec["max_ctrl_size"]
                   <= 0.10)
        ok_hi = np.isfinite(rec["pow_2"]) and rec["pow_2"] >= 0.80
        ok_lo = np.isfinite(rec["pow_half"]) and rec["pow_half"] <= 0.25
        rec["pf1_pass"] = bool(ok_size and ok_hi and ok_lo)
        pf_rows.append(rec)
    pf = pd.DataFrame(pf_rows)
    pf.to_csv(RES / "bench_pf1.csv", index=False)

    # ---------------- head-to-head calibration table --------------------
    heads = []
    methods = [("S2-bench alarm", "size_power_bench"),
               ("UCM-rho permboot", "size_power_ucm"),
               ("JS-asym permboot", "size_power_js"),
               ("Scree TW99 (S0)", "size_power_s0"),
               ("Partial F (B1)", "size_power_b1")]
    for cn in freeze["configs"]:
        for meth, col in methods:
            sel = arms[arms.config_name == cn]
            get = lambda a: (float(sel.loc[sel.arm == a, col].iloc[0])
                             if len(sel[sel.arm == a]) else np.nan)
            heads.append({
                "config_name": cn, "method": meth,
                "size_perm": get("perm_null"),
                "size_splithalf": get("splithalf"),
                "pow_half": get("pos_half"), "pow_1": get("pos_1"),
                "pow_2": get("pos_2"),
                "stat_monotone_g": None,
            })
    head_df = pd.DataFrame(heads)

    # monotonicity of continuous statistics in injected g
    mono_rows = []
    gorder = ["null", "pos_half", "pos_1", "pos_2"]
    for cn in freeze["configs"]:
        sel = arms[(arms.config_name == cn) &
                   (arms.arm.isin(gorder))].copy()
        if len(sel) < 4:
            continue
        sel["ordk"] = sel.arm.map({a: i for i, a in enumerate(gorder)})
        sel = sel.sort_values("ordk")

        def spearman(col):
            x = sel[col].to_numpy()
            rx = np.argsort(np.argsort(x))
            return float(np.corrcoef(rx, np.arange(len(x)))[0, 1])
        mono_rows.append({"config_name": cn,
                          "ucm_rho_spearman_g": spearman("ucm_rho_mean"),
                          "js_asym_spearman_g": spearman("js_asym_mean")})
    mono = pd.DataFrame(mono_rows)
    head_df = head_df.merge(mono, on="config_name", how="left")
    head_df.drop(columns=["stat_monotone_g"], inplace=True)
    head_df.to_csv(RES / "calibration_informality.csv", index=False)

    # ---------------- PF-3 / M2 ----------------------------------------
    m2 = arms[(arms.arm.isin(["m2_null", "m2_pos"]))].copy()
    m2.to_csv(RES / "bench_m2.csv", index=False)
    pf3 = {}
    try:
        pos = m2[m2.arm == "m2_pos"].iloc[0]
        pf3 = {
            "tau_ols_err": pos.get("tau_ols_err_mean"),
            "tau_trim_err": pos.get("tau_trim_onatski_err_mean"),
            "tau_trim_sign_ok": pos.get("tau_trim_onatski_sign_ok"),
            "ratio_ols_over_trim": (
                pos["tau_ols_err_mean"] / pos["tau_trim_onatski_err_mean"]
                if pos.get("tau_trim_onatski_err_mean") else np.nan),
            "pf3_pass": bool(
                pos.get("tau_trim_onatski_sign_ok", 0) >= 0.95
                and pos["tau_ols_err_mean"] >= 2.0 *
                max(pos.get("tau_trim_onatski_err_mean"), 1e-9)),
        }
    except (IndexError, KeyError):
        pass

    # ---------------- sensitivity snapshot ------------------------------
    sens = arms[arms.arm.isin(["align_top", "align_weak", "rinj_minus",
                               "rinj_plus", "hetero_eps"])].copy()
    sens.to_csv(RES / "bench_sensitivity.csv", index=False)

    # ---------------- mechanical G4 verdict -----------------------------
    main_cfgs = ["A_main", "B_main", "C_main"]
    pf1_main = pf[pf.config_name.isin(main_cfgs)]
    families_pass = int(pf1_main["pf1_pass"].sum())
    pf2_all = bool(df[df.arm == "perm_null"]["rej_s0_tw99"].mean() >= 0.98)
    verdict = {
        "G4_question": "credible result on real geometry",
        "families_pf1_pass": families_pass,
        "pf1_by_config": {r.config_name: bool(r.pf1_pass)
                          for r in pf.itertuples()},
        "pf1_primary_note":
            "PF-1 requires the straddle on the MAIN config of each family",
        "pf2_scree_false_alarm_unmodified_and_null": pf2_all,
        "pf3": pf3,
        "go_condition": {
            "null_rejection_in_[0.02,0.10]_on_geq2_families":
                bool(families_pass >= 2),
            "positive_power_straddles_predicted_frontier":
                bool(families_pass >= 2),
            "one_predeclared_finding_survives": None,
        },
    }
    surviving = []
    if families_pass >= 2:
        surviving.append("PF-1")
    if pf2_all:
        surviving.append("PF-2")
    if pf3.get("pf3_pass"):
        surviving.append("PF-3")
    verdict["go_condition"]["one_predeclared_finding_survives"] = \
        bool(surviving)
    verdict["surviving_findings"] = surviving
    verdict["verdict"] = "GO" if (families_pass >= 2 and surviving) else \
        "REVIEW"
    (RES / "bench_g4_verdict.json").write_text(json.dumps(verdict,
                                                          indent=1))

    # ---------------- figure -------------------------------------------
    make_figure(freeze, arms, main_cfgs)
    print(json.dumps(verdict, indent=1))


def make_figure(freeze, arms, main_cfgs):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    for ax, cn in zip(axes, main_cfgs):
        gs = freeze["configs"][cn]["g_star"]
        sel = arms[arms.config_name == cn].set_index("arm")
        xs, ys = [], []
        for tag, mult in [("pos_half", 0.5), ("pos_1", 1.0), ("pos_2", 2.0)]:
            if tag in sel.index:
                xs.append(mult)
                ys.append(float(sel.loc[tag, "size_power_bench"]))
        ax.plot([0] + xs, [float(sel.loc["perm_null", "size_power_bench"])]
                + ys, "o-", label="S2-bench empirical")
        ax.axhline(0.05, color="gray", lw=0.8, ls=":")
        ax.axhline(0.8, color="tab:red", lw=0.8, ls="--")
        ax.axhline(0.25, color="tab:red", lw=0.8, ls=":")
        if np.isfinite(gs):
            ax.axvspan(0.5, 2.0, color="tab:blue", alpha=0.06)
            ax.axvline(gs / gs, color="k", lw=0.8)
        ax.set_title(f"{cn} (g*={gs:.2f})")
        ax.set_xlabel("injected link strength (multiples of predicted g*)")
        ax.set_ylabel("rejection rate")
        ax.set_ylim(-0.02, 1.05)
        for meth, col, mk in [("UCM", "size_power_ucm", "s"),
                              ("JS-asym", "size_power_js", "^"),
                              ("Scree", "size_power_s0", "v")]:
            yy = [float(sel.loc[t, col]) if t in sel.index else np.nan
                  for t in ["pos_half", "pos_1", "pos_2"]]
            ax.plot(xs, yy, mk + "--", ms=4, alpha=0.6, label=meth)
        ax.legend(fontsize=7, loc="upper left")
    fig.suptitle("Benchmark frontier check: calibrated alarm vs mandatory "
                 "baselines (dashed red = PF-1 window)")
    fig.tight_layout()
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / "benchmark_frontier_check.pdf")
    print("wrote figure", FIG / "benchmark_frontier_check.pdf")


if __name__ == "__main__":
    main()
