"""SCF Phase 3: ingest Colab benchmark shards (scf_bench_<config>.zip).

Usage: python3 code/consolidate_bench.py <folder-with-zips>

Verifies each archive's manifest.json sha256 entries, then copies payloads
into data/benchmarks/raw/<config>/<arm>.parquet (+ _means.npz). Never
overwrites a LARGER existing file with a smaller one (lost-update guard,
deviation D-B1 policy). Prints a completeness table afterwards.
"""
from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "benchmarks" / "raw"


def main():
    folder = Path(sys.argv[1])
    zips = sorted(folder.glob("scf_bench_*.zip"))
    if not zips:
        print("no scf_bench_*.zip found in", folder)
        return 1
    for zp in zips:
        import zipfile

        with zipfile.ZipFile(zp) as zf:
            names = zf.namelist()
            if "manifest.json" not in names:
                print(f"[skip] {zp.name}: no manifest")
                continue
            manifest = json.loads(zf.read("manifest.json"))
            bad = []
            for arc, want in manifest.items():
                data = zf.read(arc)
                got = hashlib.sha256(data).hexdigest()
                if got != want:
                    bad.append(arc)
            if bad:
                print(f"[abort] {zp.name}: checksum mismatch {bad[:3]}")
                continue
            n_in = 0
            for arc, want in manifest.items():
                dst = RAW.parent / arc.replace("state/", "raw/")
                # archives store members relative to '.', i.e. state/<cfg>/...
                dst = ROOT / "data" / "benchmarks" / arc
                dst.parent.mkdir(parents=True, exist_ok=True)
                data = zf.read(arc)
                if dst.exists() and dst.stat().st_size > len(data):
                    continue  # keep larger local version
                tmp = dst.with_suffix(dst.suffix + ".ingest.tmp")
                tmp.write_bytes(data)
                tmp.replace(dst)
                n_in += 1
            print(f"[ok] {zp.name}: {n_in} files ingested")
    # completeness table
    import pandas as pd

    targets = {"null": 600, "perm_null": 400, "pos_half": 400, "pos_1": 400,
               "pos_2": 400, "splithalf": 300, "align_top": 200,
               "align_weak": 200, "rinj_minus": 200, "rinj_plus": 200,
               "hetero_eps": 200, "m2_null": 300, "m2_pos": 300}
    print("\ncompleteness:")
    total_short = 0
    for cfg_dir in sorted(RAW.glob("*/")):
        cfg = cfg_dir.name
        for f in sorted(cfg_dir.glob("*.parquet")):
            arm = f.stem
            t = targets.get(arm)
            if t is None:
                continue
            try:
                have = pd.read_parquet(f, columns=["rep"])["rep"].nunique()
            except Exception:
                have = -1
            mark = "OK" if have >= t else f"SHORT ({have}/{t})"
            if have < t:
                total_short += 1
            print(f"  {cfg}/{arm}: {mark}")
    print(f"\nshort cells: {total_short}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
