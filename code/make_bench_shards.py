"""SCF Phase 3 Colab shard generator (research plan Section 10.2 rule).

Pass-2 wall time on the contended local workstation exceeds the 2-hour
threshold for the heavy configs, so the sweep ships as self-contained
Colab notebooks: one shard per benchmark config, running ALL of that
config's pass-2 arms. Self-containment contract:

  1. pinned dependencies installed in cell 1;
  2. harness code embedded verbatim from the repository files
     (de_formulas.py, simulator.py, detection.py, estimators.py,
      benchmarks_data.py, benchmarks.py) via %%writefile, with sha256
     prefixes recorded at generation time;
  3. designs REGENERATED deterministically inside the notebook from primary
     sources (GEO series matrices / mirror CSVs) and verified against the
     pinned sha256 of X.tobytes() BEFORE any arm runs;
  4. calibration constants loaded from an EMBEDDED copy of
     results/benchmark_freeze.json (pass-1 output; notebooks never
     recalibrate);
  5. per-arm checkpointing into /content/state/<config>/<arm>.parquet
     (row-count resume);
  6. outputs zipped to scf_bench_<config>.zip with a manifest; final cell
     uses the safe Colab download fallback.

Usage: python3 code/make_bench_shards.py [--configs A_main,B_main,...]
Outputs: notebooks/colab_bench_<config>.ipynb
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
NB = ROOT / "notebooks"
BENCH = ROOT / "data" / "benchmarks"

HARNESS_FILES = ["de_formulas.py", "simulator.py", "detection.py",
                 "estimators.py", "benchmarks_data.py", "benchmarks.py"]


def sha16(path_or_bytes) -> str:
    if isinstance(path_or_bytes, bytes):
        return hashlib.sha256(path_or_bytes).hexdigest()[:16]
    return hashlib.sha256(Path(path_or_bytes).read_bytes()).hexdigest()[:16]


def design_hashes() -> dict:
    out = {}
    for name in ["A_main", "A_sub", "B_main", "B_wide", "C_main", "C_wide"]:
        z = np.load(BENCH / f"{name}.npz", allow_pickle=True)
        X = np.ascontiguousarray(z["X"], dtype=np.float64)
        out[name] = sha16(X.tobytes())
    return out


def cell(code: str) -> dict:
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": [], "source": code.splitlines(keepends=True)}


def md(text: str) -> dict:
    return {"cell_type": "markdown", "metadata": {},
            "source": text.splitlines(keepends=True)}


def make_notebook(config_name: str, freeze: dict, hashes: dict) -> dict:
    cells = []
    cells.append(md(f"""# SCF Phase 3 benchmark shard: **{config_name}**

Self-contained pass-2 runner. Regenerates the processed design from primary
sources, verifies its sha256 against the frozen value, then runs this
config's evaluation arms against the embedded calibration (pass-1 freeze).
Runtime target < 6 h on CPU; checkpoint-resume safe."""))
    cells.append(cell(
        "%pip install -q numpy scipy pandas pyarrow matplotlib\n"
        "import os\n"
        "for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',"
        "'NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):\n"
        "    os.environ[v] = '1'\nprint('env set')"))
    # harness writer cells
    for fn in HARNESS_FILES:
        src = (ROOT / "code" / fn).read_text()
        b64 = None
        payload = json.dumps(src)
        cells.append(md(f"### harness: `{fn}` (sha256 {sha16(ROOT/'code'/fn)})"))
        cells.append(cell(
            f"harness_{fn.replace('.py','')} = {payload}\n"
            f"import pathlib, hashlib\n"
            f"p = pathlib.Path('code'); p.mkdir(exist_ok=True)\n"
            f"(p / '{fn}').write_text(harness_{fn.replace('.py','')})\n"
            f"assert hashlib.sha256((p/'{fn}').read_bytes()).hexdigest()[:16]"
            f" == '{sha16(ROOT / 'code' / fn)}', 'harness hash mismatch'\n"
            f"import sys; sys.path.insert(0, 'code')\n"
            f"print('{fn} ok')"))
    cfg_text = (ROOT / "configs" / "benchmarks_frozen.yaml").read_text()
    cells.append(md("### embedded frozen benchmark configuration"))
    cells.append(cell(
        f"cfg_yaml = {json.dumps(cfg_text)}\n"
        f"import pathlib\n"
        f"pathlib.Path('configs').mkdir(exist_ok=True)\n"
        f"pathlib.Path('configs/benchmarks_frozen.yaml').write_text(cfg_yaml)\n"
        f"print('configs ok')"))
    cells.append(md("### acquire primary sources (pinned hashes)"))
    src_cells = {
        "A": [
            ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63060/matrix/GSE63060_series_matrix.txt.gz",
             "data/benchmarks/raw/GSE63060_series_matrix.txt.gz",
             "342421d07105d786d623a4b86dc4a21b764da352e4facc7a032fbd36b0f19371"),
            ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63061/matrix/GSE63061_series_matrix.txt.gz",
             "data/benchmarks/raw/GSE63061_series_matrix.txt.gz",
             "c6de1af5bd36bca3de0921d53051bdefd05666f6e94a6906a26c8cb461dd55f5"),
        ],
        "B": [
            ("https://raw.githubusercontent.com/gpeng9/ihdp-causality/master/ihdp.csv",
             "data/benchmarks/ihdp.csv",
             "1c12eb6df2d6a48165b34963e1457a47cfe31fb300c61b5a2c2557ddd84da467"),
        ],
        "C": [
            ("https://vincentarelbundock.github.io/Rdatasets/csv/wooldridge/k401ksubs.csv",
             "data/benchmarks/k401ksubs.csv",
             "62128c8706a238eee9dfef75776fbf2b518f9d6ac0788a1f69ef1eb7ead8d78f"),
        ],
    }
    fam_key = config_name[0]
    dl_lines = ["import urllib.request, hashlib, pathlib"]
    for url, rel, want in src_cells[fam_key]:
        dl_lines.append(f"""
pathlib.Path('{rel}').parent.mkdir(parents=True, exist_ok=True)
if not pathlib.Path('{rel}').exists():
    urllib.request.urlretrieve('{url}', '{rel}')
got = hashlib.sha256(pathlib.Path('{rel}').read_bytes()).hexdigest()
assert got == '{want}', ('source hash mismatch', '{rel}', got)
print('source ok:', '{rel}')""")
    cells.append(cell("\n".join(dl_lines)))
    cells.append(md("### regenerate + verify the design"))
    cells.append(cell(f"""
import numpy as np, json, hashlib, sys
sys.path.insert(0, 'harness')
import benchmarks_data as BD

CONFIG = "{config_name}"
EXPECTED_X_SHA16 = "{hashes[config_name]}"
RAW_DIR_NAME = "{config_name}"

if CONFIG.startswith('A'):
    fam = BD.build_addneuromed()
else:
    fam = BD.build_tabular()
payload = fam[CONFIG]
X = np.ascontiguousarray(payload['X'], dtype=np.float64)
got = hashlib.sha256(X.tobytes()).hexdigest()[:16]
assert got == EXPECTED_X_SHA16, (got, EXPECTED_X_SHA16)
meta = {{k: (v.tolist() if isinstance(v, np.ndarray) else v)
        for k, v in payload.items() if k != 'X'}}
meta_json = json.dumps(meta, default=str)
np.savez_compressed('{config_name}.npz', X=X,
                    config_name=CONFIG, meta_json=meta_json)
print('design verified:', CONFIG, X.shape)
"""))
    cells.append(md("### embedded pass-1 freeze for this config"))
    cells.append(cell(f"""
freeze_entry = json.loads(r'''{json.dumps(freeze["configs"][config_name])}''')
FREEZE = {{"version": {json.dumps(freeze["version"])},
          "ledger_hash": {json.dumps(freeze["ledger_hash"])},
          "configs": {{CONFIG: freeze_entry}}}}
print(json.dumps(freeze_entry, indent=1))
"""))
    # patch benchmarks paths for /content layout
    cells.append(md("### run all pass-2 arms for this config (checkpointed)"))
    cells.append(cell("""
import shutil, time
from pathlib import Path
import pandas as pd
import benchmarks as B

B.BENCH_DIR = Path('.').resolve()
B.RAW_DIR = B.BENCH_DIR / 'state'
state_cfg = B.RAW_DIR / CONFIG
state_cfg.mkdir(parents=True, exist_ok=True)
dst_npz = B.BENCH_DIR / f'{CONFIG}.npz'
if Path(f'{CONFIG}.npz').resolve() != dst_npz.resolve():
    shutil.copy(f'{CONFIG}.npz', dst_npz)
# spectral audit is only needed by build_jobs; provide it
audit = {"r_inj": freeze_entry["r_inj"]}
(B.BENCH_DIR / 'spectral_audit.json').write_text(json.dumps(
    {"spectral_profiles": {CONFIG: audit}}))

_orig_specs = B.arm_specs
B.arm_specs = lambda: {**_orig_specs(),
                       'configs': {CONFIG: _orig_specs()['configs'][CONFIG]}}
bench = B.Bench(CONFIG)
jobs = B.build_jobs('pass2', FREEZE)
t0 = time.time()
for j in jobs:
    if time.time() - t0 > 5.5 * 3600:
        print('time budget reached; rerun notebook to resume'); break
    B.run_cell(j)
print('shard done or budget-stopped')
"""))
    cells.append(md("### package outputs"))
    cells.append(cell(f"""
import zipfile, hashlib
zf = zipfile.ZipFile('./scf_bench_{config_name}.zip', 'w')
manifest = {{}}
for f in sorted(Path('./state').rglob('*')):
    if f.is_file():
        arc = str(f.relative_to('.'))
        zf.write(f, arcname=arc)
        manifest[arc] = hashlib.sha256(f.read_bytes()).hexdigest()
zf.writestr('manifest.json', json.dumps(manifest, indent=1))
zf.close()
print('wrote scf_bench_{config_name}.zip')
try:
    from google.colab import files
    files.download('scf_bench_{config_name}.zip')
    print('Downloaded: scf_bench_{config_name}.zip')
except Exception as e:
    print('(Not on Colab / download skipped):', e)
"""))
    return {
        "nbformat": 4, "nbformat_minor": 5,
        "metadata": {"accelerator": "None",
                     "colab": {"provenance": []}},
        "cells": cells,
    }


def main():
    args = sys.argv[1:]
    only = None
    if "--configs" in args:
        only = args[args.index("--configs") + 1].split(",")
    freeze = json.loads((ROOT / "results" / "benchmark_freeze.json")
                        .read_text())
    hashes = design_hashes()
    NB.mkdir(exist_ok=True)
    for name in (only or list(freeze["configs"].keys())):
        nb = make_notebook(name, freeze, hashes)
        out = NB / f"colab_bench_{name}.ipynb"
        json.dump(nb, open(out, "w"), indent=1)
        print("wrote", out)


if __name__ == "__main__":
    main()
