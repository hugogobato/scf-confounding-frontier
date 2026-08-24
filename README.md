# SCF: Spectral Confounding Frontier

Research code for the project "Phase Transitions for Hidden Confounding
(Detectability Frontiers and Optimal Spectral Deconfounding)".

Phase status: Phase 2 (simulation engine, gate G3) launched; grids and
thresholds frozen before any sweep data exists (see
`docs/phase2_preregistration.md`, deviation register D1-D7, and
`docs/phase2_execution_memo.md`).

## Layout

    code/     de_formulas.py (deterministic-equivalent formula sheet),
              simulator.py (DGP surface), estimators.py (frozen roster),
              detection.py (statistics + Le Cam probe), runners.py
              (9 sweep modes, checkpoint/resume), build_grids.py,
              make_shards.py (Colab notebook generator), consolidate_shards.py,
              gate_checks.py, cost_pilot.py
    configs/  frozen grid definitions grid_*.json (hash a920e554e692)
    docs/     model card, assumption ledger, preregistration, gate memos,
              detection spec with pre-data errata, novelty memo, pilot memo
    tests/    identity + behavior tests (pytest)
    notebooks/ generated self-contained Colab shards (fetch code from this
               repo at the pinned tag; jobs embedded as JSON)

## Running the Colab shards

Each notebook clones this repository at a pinned tag, executes its embedded
job list (checkpoint/resume per cell), verifies sha256 checksums, archives
outputs, and triggers the browser download. After downloading the archives:

    python3 code/consolidate_shards.py <folder-with-zips>
    python3 code/gate_checks.py

## Reproducing locally

    pip install -r requirements.txt
    python3 -m pytest tests/ -q          # 35 identity/behavior tests
    python3 code/build_grids.py          # regenerate frozen grids
    python3 code/make_shards.py          # regenerate Colab shards

## Provenance conventions

Every simulation row carries the ledger hash of the model card +
assumption ledger pair. Figures must regenerate exclusively from parquet.
Thresholds in memos were fixed before data generation; post-data changes
require a deviation entry.

License: TBD.
