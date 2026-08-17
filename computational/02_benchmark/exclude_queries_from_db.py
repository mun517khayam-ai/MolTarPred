#!/usr/bin/env python3
"""
Remove the master query ligands from a version's reference database, so no version
gets a self-match advantage another lacks. Run this on EVERY version's grouped
database before predicting, using the SAME master query file for all of them.

Usage:
  python exclude_queries_from_db.py --db <grouped_db.csv> [--query query_master_FDA.csv] [--out <db_masked.csv>]

Default --out is "<db stem>_masked.csv" next to the input.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

CODE = "/Users/macbook/chembl/Code"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True, help="grouped DB csv (Ligand-id, SMILES, Target-id)")
    ap.add_argument("--query", default=f"{CODE}/query_master_FDA.csv")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    db = pd.read_csv(args.db, dtype=str)
    q = pd.read_csv(args.query, dtype=str)
    qset = set(q["Ligand-id"].astype(str))

    before = len(db)
    masked = db[~db["Ligand-id"].astype(str).isin(qset)]
    removed = before - len(masked)

    out = Path(args.out) if args.out else Path(args.db).with_name(Path(args.db).stem + "_masked.csv")
    masked.to_csv(out, index=False)
    print(f"[done] {Path(args.db).name}: removed {removed} query ligands "
          f"({before} -> {len(masked)}); wrote {out}")


if __name__ == "__main__":
    main()
