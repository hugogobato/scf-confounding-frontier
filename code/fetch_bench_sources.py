"""SCF Phase 3: fetch primary benchmark sources with pinned checksums.

Usage: python3 code/fetch_bench_sources.py --family A|B|C [--all]

Downloads only what the requested family needs into data/benchmarks/
(matching benchmarks_data.py's expected paths) and verifies sha256 BEFORE
any preprocessing runs. Sources already on disk and hash-valid are kept.
"""
from __future__ import annotations

import argparse
import hashlib
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "benchmarks"

SOURCES = {
    "A": [(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63060/matrix/GSE63060_series_matrix.txt.gz",
        "raw/GSE63060_series_matrix.txt.gz",
        "342421d07105d786d623a4b86dc4a21b764da352e4facc7a032fbd36b0f19371"),
        ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63061/matrix/GSE63061_series_matrix.txt.gz",
         "raw/GSE63061_series_matrix.txt.gz",
         "c6de1af5bd36bca3de0921d53051bdefd05666f6e94a6906a26c8cb461dd55f5")],
    "B": [(
        "https://raw.githubusercontent.com/gpeng9/ihdp-causality/master/ihdp.csv",
        "ihdp.csv",
        "1c12eb6df2d6a48165b34963e1457a47cfe31fb300c61b5a2c2557ddd84da467")],
    "C": [(
        "https://vincentarelbundock.github.io/Rdatasets/csv/wooldridge/k401ksubs.csv",
        "k401ksubs.csv",
        "62128c8706a238eee9dfef75776fbf2b518f9d6ac0788a1f69ef1eb7ead8d78f")],
}


def fetch(family: str) -> None:
    for url, rel, want in SOURCES[family]:
        dst = OUT / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.exists():
            got = hashlib.sha256(dst.read_bytes()).hexdigest()
            if got == want:
                print(f"[keep] {rel} (hash ok)")
                continue
            print(f"[refetch] {rel} (hash mismatch)")
            dst.unlink()
        tmp = dst.with_suffix(dst.suffix + ".dl")
        print(f"[down] {url}")
        urllib.request.urlretrieve(url, tmp)
        got = hashlib.sha256(tmp.read_bytes()).hexdigest()
        if got != want:
            tmp.unlink()
            raise SystemExit(f"source hash mismatch for {rel}: {got}")
        tmp.replace(dst)
        print(f"[ok] {rel}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", choices=["A", "B", "C"])
    ap.add_argument("--all", action="store_true")
    args = ap.parse_args()
    fams = ["A", "B", "C"] if args.all else [args.family]
    for f in fams:
        fetch(f)


if __name__ == "__main__":
    main()
