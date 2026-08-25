"""SCF Phase 3 Colab shard driver (clone-at-tag mode).

Usage inside a cloned repository pinned at phase3-freeze:

    python3 code/colab_shard_main.py CONFIG [--time-budget SECONDS]

Steps:
  1. fetch the primary sources for CONFIG's family (pinned sha256),
  2. regenerate the processed design deterministically and verify its
     sha256 prefix against configs/benchmarks_design_hashes.json,
  3. load the frozen calibration (results/benchmark_freeze.json),
  4. run ALL pass-2 arms for this config with checkpoint-resume into
     ./state/<config>/<arm>.parquet,
  5. package outputs into scf_bench_<config>.zip with a sha256 manifest.

The notebook then only needs to offer the download fallback.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    import os

    os.environ.setdefault(_v, "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--time-budget", type=float, default=5.5 * 3600)
    args = ap.parse_args()
    config = args.config

    import numpy as np

    import benchmarks_data as BD
    import benchmarks as B

    # ---- 1. primary sources ------------------------------------------
    fam_key = config[0]
    sys.argv = ["fetch", "--family", fam_key]
    try:
        BD_fetch = __import__("fetch_bench_sources")
        BD_fetch.main()
    except SystemExit:
        pass

    # ---- 2. regenerate + verify design -------------------------------
    hashes = json.loads(
        (ROOT / "configs" / "benchmarks_design_hashes.json").read_text())
    expected = hashes[config]
    if config.startswith("A"):
        fam = BD.build_addneuromed()
    else:
        fam = BD.build_tabular()
    payload = fam[config]
    X = np.ascontiguousarray(payload["X"], dtype=np.float64)
    got = hashlib.sha256(X.tobytes()).hexdigest()[:16]
    if got != expected:
        raise SystemExit(f"design hash mismatch for {config}: {got}")
    meta = {k: (v.tolist() if isinstance(v, np.ndarray) else v)
            for k, v in payload.items() if k != "X"}
    B.BENCH_DIR = ROOT / "data" / "benchmarks"
    B.BENCH_DIR.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(B.BENCH_DIR / f"{config}.npz", X=X,
                        config_name=config,
                        meta_json=json.dumps(meta, default=str))
    print(f"[design] {config} verified {X.shape}")

    # ---- 3. frozen calibration ---------------------------------------
    freeze = json.loads(
        (ROOT / "results" / "benchmark_freeze.json").read_text())
    entry = freeze["configs"][config]
    coord_scales = np.asarray(entry["coord_scales"], float)

    # ---- 4. run all pass-2 arms for this config ----------------------
    B.RAW_DIR = ROOT / "state"
    state_cfg = B.RAW_DIR / config
    state_cfg.mkdir(parents=True, exist_ok=True)
    (B.BENCH_DIR / "spectral_audit.json").write_text(json.dumps(
        {"spectral_profiles": {config: {"r_inj": entry["r_inj"]}}}))

    bench = B.Bench(config)
    _orig_specs = B.arm_specs
    B.arm_specs = lambda: {**_orig_specs(),
                           "configs": {config:
                                       _orig_specs()["configs"][config]}}
    jobs = B.build_jobs("pass2", freeze)
    t0 = time.time()
    interrupted = False
    for j in jobs:
        if time.time() - t0 > args.time_budget:
            print("[budget] stopping; rerun to resume", flush=True)
            interrupted = True
            break
        B.run_cell(j)

    # ---- 5. package ---------------------------------------------------
    zf_path = ROOT / f"scf_bench_{config}.zip"
    import zipfile

    with zipfile.ZipFile(zf_path, "w") as zf:
        manifest = {}
        for f in sorted((ROOT / "state").rglob("*")):
            if not f.is_file():
                continue
            arc = str(f.relative_to(ROOT))
            zf.write(f, arcname=arc)
            manifest[arc] = hashlib.sha256(f.read_bytes()).hexdigest()
        zf.writestr("manifest.json", json.dumps(manifest, indent=1))
    print(f"[done] wrote {zf_path.name}"
          + (" (INCOMPLETE - budget)" if interrupted else " (complete)"),
          flush=True)


if __name__ == "__main__":
    main()
