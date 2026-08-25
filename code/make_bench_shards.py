"""SCF Phase 3 Colab shard generator (clone-at-tag mode, v2).

Fixes vs v1 (2026-08-25):
  * REMOVED the invalid metadata accelerator key ("Unknown accelerator:
    None" Colab error) - CPU runtime is Colab's default when the key is
    absent.
  * Notebooks no longer embed harness code: they CLONE this repository at
    the pinned tag (phase2 pattern) and verify harness sha256 prefixes
    against generation-time values before running anything.
  * All logic lives in code/colab_shard_main.py inside the repo; the
    notebook is three functional cells (install, clone+verify, run+download).

Usage: python3 code/make_bench_shards.py [--tag phase3-freeze] [--configs ...]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NB = ROOT / "notebooks"

HARNESS_FILES = ["code/de_formulas.py", "code/simulator.py",
                 "code/detection.py", "code/estimators.py",
                 "code/benchmarks_data.py", "code/fetch_bench_sources.py",
                 "code/benchmarks.py", "code/colab_shard_main.py",
                 "code/bench_freeze.py", "code/build_grids.py"]


def sha16(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16]


def cell(code: str) -> dict:
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": [], "source": code.splitlines(keepends=True)}


def md(text: str) -> dict:
    return {"cell_type": "markdown", "metadata": {},
            "source": text.splitlines(keepends=True)}


def make_notebook(config: str, tag: str) -> dict:
    checks = "\n".join(
        f"assert sha('{f}') == '{sha16(ROOT / f)}', 'harness mismatch: {f}'"
        for f in HARNESS_FILES)
    cells = [
        md(f"""# SCF Phase 3 benchmark shard: **{config}**

Clone-at-tag runner (pinned `{tag}`). Runtime: CPU (default). Checkpoint-
resume safe: if the session dies, Runtime > Run all resumes from saved
state. Total wall target < 6 h."""),
        cell("%pip install -q numpy scipy pandas pyarrow\n"
             "print('deps ok')"),
        cell(f"""import os, subprocess, hashlib
from pathlib import Path
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
          'NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[v] = '1'
TAG = "{tag}"
REPO = "https://github.com/hugogobato/scf-confounding-frontier.git"
if not Path("scf").exists():
    subprocess.run(["git", "clone", "--depth", "1", "-b", TAG, REPO, "scf"],
                   check=True)
os.chdir("scf")
def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
{checks}
print("harness verified at tag", TAG)"""),
        cell(f"""CONFIG = "{config}"
!python code/colab_shard_main.py {{CONFIG}} --time-budget 18000"""),
        cell(f"""try:
    from google.colab import files
    files.download("scf_bench_{config}.zip")
    print("Downloaded: scf_bench_{config}.zip")
except Exception as e:
    print("(Not on Colab / download skipped):", e)
print("If incomplete: Runtime > Run all to resume from state/")"""),
    ]
    return {"nbformat": 4, "nbformat_minor": 5,
            "metadata": {"colab": {"provenance": []},
                         "kernelspec": {"name": "python3",
                                        "display_name": "Python 3"},
                         "language_info": {"name": "python"}},
            "cells": cells}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="phase3-freeze")
    ap.add_argument("--configs")
    args = ap.parse_args()
    # refuse to pin a tag that does not exist or does not point at HEAD
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                          capture_output=True, text=True).stdout.strip()
    tgt = subprocess.run(["git", "rev-parse", args.tag + "^{{commit}}"],
                         cwd=ROOT, capture_output=True, text=True)
    if tgt.returncode != 0:
        raise SystemExit(f"tag {args.tag} does not exist - create it first "
                         f"(git tag {args.tag} && git push origin {args.tag})")
    if tgt.stdout.strip() != head:
        raise SystemExit(f"tag {args.tag} points at "
                         f"{tgt.stdout.strip()[:12]}, not HEAD "
                         f"{head[:12]} - move the tag first")
    configs = (args.configs.split(",") if args.configs else
               [p.stem.replace("colab_bench_", "")
                for p in sorted(NB.glob("colab_bench_*.ipynb"))])
    NB.mkdir(exist_ok=True)
    for config in configs:
        nb = make_notebook(config, args.tag)
        out = NB / f"colab_bench_{config}.ipynb"
        json.dump(nb, open(out, "w"), indent=1)
        print("wrote", out)


if __name__ == "__main__":
    main()
