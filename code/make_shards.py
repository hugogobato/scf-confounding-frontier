"""Phase 2 Colab shard builder.

Packs frozen grid jobs (configs/grid_*.json) into <= 36 self-contained
notebooks targeting <= 5.0 projected hours each (hard stop logic inside),
using a cost model calibrated on WP 1.5 timings (same workstation class):

    t_rep(s) ~ OVERHEAD + ALPHA * n * p * min(n, p) + BETA * min(n, p)^3

with ALPHA = 1.56e-10, BETA = 6.9e-10, OVERHEAD = 0.22 fitted to the pilot
log points (2000,400)->0.27s, (2000,2000)->~7s, (2000,10000)->~12s, and
mode/safety multipliers below. All packing happens BEFORE any data is
generated; cells are atomic (never split across shards) and resume-safe.

Outputs: notebooks/colab_shard_XX.ipynb + data/sim/shard_manifest.json.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
NB_DIR = ROOT / "notebooks"
SIM = ROOT / "data" / "sim"

ALPHA, BETA, OVERHEAD = 1.56e-10, 6.9e-10, 0.22
MODE_MULT = {
    "correctness": 1.15, "nullcal": 1.0, "power": 1.0, "alignment": 1.0,
    "estimation": 3.0, "crossover": 2.2, "robustness": 2.2, "m2": 1.2,
    "scaling": 1.15,
}
TWIN_MULT = 1.8          # gamma=0 twin arms refit the roster (c > 1 cells)
SAFETY = 1.8             # Colab vCPU vs loaded-local calibration
WORKERS = 2              # Colab notebook workers; wall = core / WORKERS
TARGET_H, HARD_H = 5.0, 5.5
MAX_SHARDS = 40
SWEEP_ORDER = ["correctness", "nullcal", "estimation", "crossover", "power",
               "alignment", "robustness", "m2", "scaling"]
CODE_FILES = ["de_formulas.py", "simulator.py", "estimators.py",
              "detection.py", "runners.py"]
# Code is fetched from GitHub at a pinned tag (public repo), not embedded;
# regenerate notebooks + retag whenever code changes.
GITHUB_REPO = "hugogobato/scf-confounding-frontier"
PINNED_TAG = "phase2-freeze"


def rep_cost(cfg: dict, reps: int, mode: str) -> float:
    """Core-seconds for a cell."""
    n, p = cfg["n"], cfg["p"]
    m = min(n, p)
    base = OVERHEAD + ALPHA * n * p * m + BETA * m ** 3
    mult = MODE_MULT[mode] * (TWIN_MULT if cfg.get("twin_gamma0") else 1.0)
    return base * mult * SAFETY * reps


def shard_wall_hours(jobs: list[dict]) -> float:
    """Wall-hours at WORKERS-fold process parallelism."""
    return sum(rep_cost(j["config"], j["reps"], j["mode"])
               for j in jobs) / WORKERS / 3600.0


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_notebook(shard_id: int, jobs: list[dict], sweeps: list[str],
                   proj_h: float, grids_hash: str, code_hashes: dict) -> str:
    nb = {
        "cells": [],
        "metadata": {"colab": {"provenance": []},
                     "kernelspec": {"name": "python3",
                                    "display_name": "Python 3"}},
        "nbformat": 4, "nbformat_minor": 0,
    }

    def md(text):
        nb["cells"].append({"cell_type": "markdown", "metadata": {},
                            "source": text.splitlines(keepends=True)})

    def code(src):
        nb["cells"].append({"cell_type": "code", "metadata": {},
                            "execution_count": None, "outputs": [],
                            "source": src.splitlines(keepends=True)})

    md(f"""# SCF Phase 2 shard {shard_id:02d}

Sweeps: {', '.join(sweeps)}. Jobs: {len(jobs)}. Projected: {proj_h:.1f} h
(safety-factor {SAFETY}x applied). Grids hash: `{grids_hash}`.
Code source: github.com/{GITHUB_REPO} @ tag `{PINNED_TAG}` (pinned for
reproducibility).
Pre-registration: `docs/phase2_preregistration.md` (thresholds frozen before
any data generation; deviation register D1-D7 included there).

Resume-safe: completed cells are skipped on rerun (checkpoint parquet per
cell). If the notebook approaches the Colab wall limit it finishes the
current cell and stops cleanly; rerun to continue.""")
    code("""import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"
!pip install -q "numpy>=2.0" "scipy>=1.14" "pandas>=2.2" "pyarrow>=16" scikit-learn""")
    code(f"""!git clone --depth 1 --branch {PINNED_TAG} \\
    https://github.com/{GITHUB_REPO}.git scf_repo
import sys, hashlib, json
sys.path.insert(0, "scf_repo/code")
# verify the pinned code matches the manifest recorded at generation time
EXPECTED = json.loads({json.dumps(json.dumps(code_hashes))})
for fname, short in EXPECTED.items():
    h = hashlib.sha256(open(f"scf_repo/code/{{fname}}", "rb").read()).hexdigest()[:12]
    assert h == short, f"code mismatch: {{fname}} ({{h}} != {{short}})"
print("code verified against generation-time hashes")""")
    code(
        "import json, time, traceback\n"
        "from multiprocessing import Pool\n"
        "from runners import run_cell\n\n"
        f"JOBS = json.loads({json.dumps(json.dumps(jobs))})\n\n"
        "def _safe(job):\n"
        "    try:\n"
        "        return run_cell(job)\n"
        "    except Exception as e:\n"
        "        print('[FAIL]', job['config_id'], repr(e))\n"
        "        traceback.print_exc()\n"
        "        return job['config_id'], -1.0\n\n"
        "t0 = time.time()\n"
        "results = []\n"
        "for i, job in enumerate(JOBS):\n"
        "    if time.time() - t0 > 8.6 * 3600:\n"
        "        print('[WALL LIMIT] stopping cleanly after', i, 'jobs')\n"
        "        break\n"
        "    results.append(_safe(job))\n"
        "print('shard done:', results)"
    )
    code(
        "import hashlib, json, glob, os\n"
        "manifest = {'shard_id': " + str(shard_id) + ", 'files': {}}\n"
        "os.makedirs('data', exist_ok=True)\n"
        "for f in sorted(glob.glob('data/**/*.parquet', recursive=True)) + \\\n"
        "         sorted(glob.glob('data/**/*.npz', recursive=True)):\n"
        "    h = hashlib.sha256(open(f, 'rb').read()).hexdigest()\n"
        "    manifest['files'][f] = h\n"
        "with open('data/manifest.json', 'w') as fh:\n"
        "    json.dump(manifest, fh, indent=1)\n"
        "print(json.dumps(manifest['files'], indent=1))"
    )
    code(
        "import shutil\n"
        "archive = shutil.make_archive('scf_shard_{:02d}'.format(" +
        str(shard_id) + "), 'zip', 'data')\n"
        "print('archived:', archive)\n"
        "output_file = archive\n"
        "try:\n"
        "    from google.colab import files\n"
        "    files.download(output_file)\n"
        "    print('Downloaded:', output_file)\n"
        "except Exception as e:\n"
        "    print('(Not on Colab / download skipped):', e)"
    )
    import json as _j

    return _j.dumps(nb, indent=1)


def main():
    grids_hash = sha256_file(ROOT / "configs" / "grid_correctness.json")[:12]
    code_hashes = {f: sha256_file(CODE / f)[:12] for f in CODE_FILES}
    jobs_by_sweep = {}
    for sweep in SWEEP_ORDER:
        gp = ROOT / "configs" / f"grid_{sweep}.json"
        jobs = json.loads(gp.read_text())
        for j in jobs:
            j["raw_path"] = j["raw_path"]
        jobs_by_sweep[sweep] = jobs

    # pack greedily: sweeps in priority order, biggest cells first; a bin may
    # mix sweeps only when absorbing leftovers (first-fit by wall time)
    all_jobs = []
    for rank, sweep in enumerate(SWEEP_ORDER):
        for j in jobs_by_sweep[sweep]:
            j["_rank"] = rank
            j["_sweep"] = sweep
            all_jobs.append(j)
    bins: list[list[dict]] = []
    bins_core: list[float] = []
    for sweep in SWEEP_ORDER:
        jobs = sorted((j for j in all_jobs if j["_sweep"] == sweep),
                      key=lambda j: -rep_cost(j["config"], j["reps"],
                                              j["mode"]))
        cur, cur_core = [], 0.0
        for j in jobs:
            s_core = rep_cost(j["config"], j["reps"], j["mode"])
            if cur and (cur_core + s_core) / WORKERS > TARGET_H * 3600:
                bins.append(cur)
                bins_core.append(cur_core)
                cur, cur_core = [], 0.0
            cur.append(j)
            cur_core += s_core
        if cur:
            # try to absorb the leftover into an existing open-size bin
            for bi, bc in enumerate(bins):
                if (bins_core[bi] + cur_core) / WORKERS <= TARGET_H * 3600 \
                        and len(bins[bi]) + len(cur) <= 40:
                    bins[bi].extend(cur)
                    bins_core[bi] += cur_core
                    cur, cur_core = [], 0.0
                    break
            if cur:
                bins.append(cur)
                bins_core.append(cur_core)
    shards = [(sorted({j["_sweep"] for j in b}), b, c)
              for b, c in zip(bins, bins_core)]
    over = [(sweeps, c / WORKERS / 3600) for sweeps, _, c in shards
            if c / WORKERS > HARD_H * 3600]
    print(f"{len(shards)} shards packed; over-hard-cap: {over}")

    NB_DIR.mkdir(exist_ok=True)
    SIM.mkdir(parents=True, exist_ok=True)
    manifest = {"grids_hash": grids_hash, "code_hashes": code_hashes,
                "github_repo": GITHUB_REPO, "pinned_tag": PINNED_TAG,
                "cost_model": {"alpha": ALPHA, "beta": BETA,
                               "overhead": OVERHEAD, "safety": SAFETY},
                "shards": []}
    for i, (sweeps, jobs, secs) in enumerate(shards, start=1):
        proj_h = secs / WORKERS / 3600
        nb_text = build_notebook(i, jobs, sweeps, proj_h, grids_hash,
                                 code_hashes)
        path = NB_DIR / f"colab_shard_{i:02d}.ipynb"
        path.write_text(nb_text)
        for j in jobs:
            j["shard"] = path.name
        manifest["shards"].append({
            "notebook": path.name, "sweeps": sweeps, "jobs": len(jobs),
            "projected_h": round(proj_h, 2),
            "config_ids": [j["config_id"] for j in jobs],
        })
        print(f"{path.name}: {len(jobs)} jobs, {proj_h:.1f} h projected")
    with open(SIM / "shard_manifest.json", "w") as f:
        json.dump(manifest, f, indent=1)
    total_h = sum(s["projected_h"] for s in manifest["shards"])
    print(f"total projected WALL: {total_h:.0f} h across "
          f"{len(manifest['shards'])} notebooks -> {SIM/'shard_manifest.json'}")


if __name__ == "__main__":
    main()
