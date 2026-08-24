"""Consolidate returned Colab shard archives into sweep-level parquets.

Usage:
    python3 code/consolidate_shards.py <shard_zip_or_dir> [<more>...]

Accepts the downloaded scf_shard_XX.zip archives (or a directory of them).
For each: verify sha256 against its embedded manifest.json, copy payloads
into data/sim/<sweep>/raw|means/, then rebuild per-sweep consolidated
parquets and completeness reports against configs/grid_*.json.
"""
from __future__ import annotations

import hashlib
import json
import shutil
import sys
import tempfile
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SIM = ROOT / "data" / "sim"


def fresh_cid(cfg_kwargs: dict) -> str:
    """Authoritative config id under the pinned harness.

    Erratum E4 / deviation D8: the config_id column frozen inside
    configs/grid_*.json was computed by a pre-freeze Config.cid; the pinned
    tag's cid function is authoritative (it also seeds every rep). Grid
    dicts are unchanged; only the id layer is reconciled here.
    """
    import sys
    sys.path.insert(0, str(ROOT / "code"))
    from simulator import Config
    return Config(**cfg_kwargs).cid


def write_remap():
    """Audit file: stale grid id -> fresh id for every grid job."""
    remap = {}
    for gp in sorted((ROOT / "configs").glob("grid_*.json")):
        for j in json.loads(gp.read_text()):
            remap[j["config_id"]] = fresh_cid(j["config"])
    out = ROOT / "configs" / "cid_remap.json"
    out.write_text(json.dumps(remap, indent=1, sort_keys=True))
    return len(remap)


def sha256_file(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _locate(td: Path, key: str):
    """Manifest keys are /content-relative ('data/sim/...'); archives may
    store members either with or without the 'data/' prefix."""
    for cand in (td / key, td / key.removeprefix("data" + "/")):
        if cand.exists():
            return cand
    return None


def _sweep_of(key: str) -> str:
    parts = Path(key).parts
    return parts[parts.index("sim") + 1] if "sim" in parts else "unknown"


def ingest_archive(archive: Path) -> dict:
    report = {"archive": archive.name, "files_ok": 0, "files_bad": [],
              "copied": 0}
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(archive) as z:
            z.extractall(td)
        td = Path(td)
        mpath = td / "data" / "manifest.json"
        if not mpath.exists():
            mpath = td / "manifest.json"
        if not mpath.exists():
            report["error"] = "missing manifest.json"
            return report
        manifest = json.loads(mpath.read_text())
        for f, h in manifest["files"].items():
            fp = _locate(td, f)
            if fp is None or sha256_file(fp) != h:
                report["files_bad"].append(f)
                continue
            report["files_ok"] += 1
            sweep = _sweep_of(f)
            dst = SIM / sweep / ("raw" if fp.suffix == ".parquet" else "means")
            dst.mkdir(parents=True, exist_ok=True)
            shutil.copy2(fp, dst / fp.name)
            report["copied"] += 1
    return report


def consolidate_sweep(sweep: str) -> dict:
    grid_path = ROOT / "configs" / f"grid_{sweep}.json"
    if not grid_path.exists():
        return {}
    grid = json.loads(grid_path.read_text())
    raw_dir = SIM / sweep / "raw"
    parts, have = [], {}
    if raw_dir.exists():
        for f in sorted(raw_dir.glob("*.parquet")):
            df = pd.read_parquet(f)
            parts.append(df)
            cid = str(df["config_id"].iloc[0])
            have[cid] = max(have.get(cid, 0), int(df["rep"].nunique()))
    status = []
    for j in grid:
        cid = fresh_cid(j["config"])
        got = have.get(cid, 0)
        status.append({
            "config_id": cid, "sweep": sweep,
            "label": j["config"].get("profile"), "reps_expected": j["reps"],
            "reps_done": min(got, j["reps"]),
            "complete": got >= j["reps"],
        })
    st = pd.DataFrame(status)
    out = SIM / sweep
    out.mkdir(parents=True, exist_ok=True)
    st.to_csv(out / "completeness.csv", index=False)
    summary = {
        "sweep": sweep,
        "cells": len(grid),
        "cells_complete": int(st["complete"].sum()),
        "reps_expected": int(st["reps_expected"].sum()),
        "reps_done": int(st["reps_done"].sum()),
    }
    if parts:
        full = pd.concat(parts, ignore_index=True)
        full.to_parquet(out / f"{sweep}_results.parquet", index=False)
        summary["rows"] = len(full)
    return summary


def main():
    args = sys.argv[1:]
    if not args:
        print(__doc__)
        sys.exit(1)
    n = write_remap()
    print(f"cid_remap.json written ({n} ids reconciled, deviation D8)")
    for a in args:
        p = Path(a)
        if p.is_dir():
            arcs = sorted(p.glob("scf_shard_*.zip"))
        else:
            arcs = [p]
        for arc in arcs:
            rep = ingest_archive(arc)
            print(rep)
    all_summary = []
    for sweep in ("correctness", "nullcal", "power", "alignment",
                  "estimation", "crossover", "robustness", "m2", "scaling"):
        s = consolidate_sweep(sweep)
        if s:
            all_summary.append(s)
            print(s)
    (SIM / "consolidation_report.json").write_text(
        json.dumps(all_summary, indent=1))


if __name__ == "__main__":
    main()
