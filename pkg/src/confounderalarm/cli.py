"""CLI: python -m confounderalarm --csv data.csv --y outcome [--treatment t]
[--exclude col1,col2] [--n-perm 400]

Prints the alarm verdict as JSON to stdout."""
from __future__ import annotations

import argparse
import json
import sys

import pandas as pd


def main(argv=None):
    ap = argparse.ArgumentParser(prog="confounderalarm")
    ap.add_argument("--csv", required=True)
    ap.add_argument("--y", required=True, help="response column")
    ap.add_argument("--treatment", default=None,
                    help="optional treatment column (enables the "
                         "trim-then-regress adjustment)")
    ap.add_argument("--exclude", default="",
                    help="comma-separated columns to drop from the design")
    ap.add_argument("--n-perm", type=int, default=400)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=20260823)
    args = ap.parse_args(argv)

    df = pd.read_csv(args.csv)
    y = df[args.y].to_numpy(dtype=float)
    drop = {args.y}
    if args.treatment:
        drop.add(args.treatment)
    if args.exclude:
        drop |= {c.strip() for c in args.exclude.split(",") if c.strip()}
    Xdf = df.drop(columns=list(drop)).select_dtypes(include=["number"])
    if Xdf.shape[1] == 0:
        print(json.dumps({"error": "no numeric design columns"}))
        return 2
    D = df[args.treatment].to_numpy(dtype=float) if args.treatment else None
    from .core import fit_alarm

    rep = fit_alarm(y, Xdf.to_numpy(dtype=float), D=D, alpha=args.alpha,
                    n_perm=args.n_perm, seed=args.seed)
    rep["design_columns"] = list(Xdf.columns)
    print(json.dumps(dict(rep), indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
