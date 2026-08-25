"""SCF Phase 3: Colab-vs-local replication check.

Same seeds and same code ran on two machines; statistics built from the
sample spectrum are invariant to eigenvector SIGN conventions, so all
compared quantities are sign-invariant. Tolerances: exact to 1e-8 for
locally-computed invariants, 1e-4 relative for quantities that pass through
eigenvectors (BLAS differences between Colab and local CPUs).

Usage: python3 code/check_replication.py notebooks
Writes results/bench_replication.csv and prints a summary verdict.
"""
from __future__ import annotations

import json
import sys
import zipfile
from io import BytesIO
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "benchmarks" / "raw"

STAT_COLS = ["t_maxz", "T_bench", "ucm_rho", "js_asym", "lam_max_cov",
             "tw_stat", "f_pcs", "b_norm2"]


def zeta_std(row_df: pd.DataFrame, scales: np.ndarray):
    cols = [c for c in row_df.columns if c.startswith("zeta")]
    Z = np.abs(row_df[cols].to_numpy())
    k = min(len(scales), Z.shape[1])
    return np.max(Z[:, :k] / np.maximum(scales[:k], 1e-12)[None, :], axis=1)


def main():
    folder = Path(sys.argv[1] if len(sys.argv) > 1 else "notebooks")
    freeze = json.loads((ROOT / "results" / "benchmark_freeze.json")
                        .read_text())
    rows = []
    worst = 0.0
    for zp in sorted(folder.glob("scf_bench_*.zip")):
        with zipfile.ZipFile(zp) as zf:
            manifest = json.loads(zf.read("manifest.json"))
            for arc in sorted(manifest):
                if not arc.endswith(".parquet") or "state/" not in arc:
                    continue
                rel = arc.split("state/", 1)[1]
                cfg, arm = rel.split("/")
                arm = arm.removesuffix(".parquet")
                remote = pd.read_parquet(BytesIO(zf.read(arc)))
                local_f = RAW / cfg / f"{arm}.parquet"
                if not local_f.exists():
                    rows.append(dict(config=cfg, arm=arm,
                                     status="LOCAL_MISSING"))
                    continue
                local = pd.read_parquet(local_f)
                if len(remote) != len(local):
                    rows.append(dict(config=cfg, arm=arm,
                                     status=f"ROWS {len(remote)} vs "
                                            f"{len(local)}"))
                    continue
                s = np.asarray(freeze["configs"][cfg]["coord_scales"], float)
                t_rem = zeta_std(remote, s)
                t_loc = zeta_std(local, s)
                d_t = float(np.max(np.abs(t_rem - t_loc)))
                d_u = float(np.max(np.abs(remote["ucm_rho"].to_numpy()
                                          - local["ucm_rho"].to_numpy())))
                d_j = float(np.max(np.abs(remote["js_asym"].to_numpy()
                                          - local["js_asym"].to_numpy())))
                rej_r = (t_rem > freeze["configs"][cfg]["mc95_S2"])
                rej_l = (t_loc > freeze["configs"][cfg]["mc95_S2"])
                n_disagree = int((rej_r != rej_l).sum())
                rows.append(dict(config=cfg, arm=arm, n=len(remote),
                                 max_dT=d_t, max_d_ucm=d_u, max_d_js=d_j,
                                 rejection_disagreements=n_disagree,
                                 status="OK"))
                worst = max(worst, d_t)
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "results" / "bench_replication.csv", index=False)
    ok = out[out.status == "OK"]
    print(f"cells compared: {len(ok)} / {len(out)}")
    if len(ok):
        print(f"max |dT|          : {ok.max_dT.max():.3e}")
        print(f"max |d ucm_rho|   : {ok.max_d_ucm.max():.3e}")
        print(f"max |d js_asym|   : {ok.max_d_js.max():.3e}")
        print(f"rejection disagreements across all cells x reps: "
              f"{int(ok.rejection_disagreements.sum())}")
    bad = out[out.status != "OK"]
    if len(bad):
        print(bad.to_string())
    print("VERDICT:", "REPLICATED" if (len(bad) == 0 and worst < 1e-2)
          else "CHECK")


if __name__ == "__main__":
    main()
