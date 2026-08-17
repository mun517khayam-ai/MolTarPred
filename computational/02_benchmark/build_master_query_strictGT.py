#!/usr/bin/env python3
"""
Build a SECOND master query file that is identical to query_master_FDA.csv except the
ground-truth ("true target") column uses the STRICT target definition instead of the
loose one. Same 100 molecules, same SMILES -> the only variable is the answer key.

Strict definition source: v2_chembl_36_filter.csv
  (<=10000 nM, IC50/Ki/EC50, relation '=', assay_type 'B', multi/complex removed) --
the strict analog of the loose Chembl_36_filter.csv used for query_master_FDA.csv.

Output: query_master_FDA_strictGT.csv  [Ligand-id, SMILES, Target-id(strict)]
"""
from __future__ import annotations
import pandas as pd

CODE = "/Users/macbook/chembl/Code"
MASTER      = f"{CODE}/query_master_FDA.csv"                 # frozen molecules + loose GT
STRICT_FILT = f"{CODE}/v2/v2_chembl_36_filter.csv"          # strict GT source
OUT         = f"{CODE}/query_master_FDA_strictGT.csv"


def main() -> None:
    master = pd.read_csv(MASTER, dtype=str)
    ligs = set(master["Ligand-id"].astype(str))

    f = pd.read_csv(STRICT_FILT, usecols=["molecule_chembl_id", "target_chembl_id"], dtype=str)
    f["lig"] = f["molecule_chembl_id"].str.replace("CHEMBL", "", regex=False)
    f["tgt"] = f["target_chembl_id"].str.replace("CHEMBL", "", regex=False)
    f = f[f["lig"].isin(ligs)]
    strict_gt = (f.groupby("lig")["tgt"]
                   .agg(lambda s: ":".join(sorted(set(s), key=int)))
                   .to_dict())

    out = master.copy()
    out["Target-id"] = out["Ligand-id"].astype(str).map(strict_gt)
    if out["Target-id"].isna().any():
        raise SystemExit("Some molecules have no strict target — unexpected.")
    out.to_csv(OUT, index=False)

    loose_n = master["Target-id"].str.split(":").map(len)
    strict_n = out["Target-id"].str.split(":").map(len)
    print(f"[done] wrote {OUT} ({len(out)} molecules)")
    print(f"  true targets/query  LOOSE: mean={loose_n.mean():.1f} median={loose_n.median():.0f} "
          f"max={loose_n.max()}")
    print(f"  true targets/query  STRICT: mean={strict_n.mean():.1f} median={strict_n.median():.0f} "
          f"max={strict_n.max()}")
    print(f"  mean shrink factor: {strict_n.mean()/loose_n.mean():.2f}x")


if __name__ == "__main__":
    main()
