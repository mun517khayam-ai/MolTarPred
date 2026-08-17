#!/usr/bin/env python3
"""
Build a FROZEN master query set for fair comparison across all MolTarPred versions
(V1 fingerprints/loose DB, V2 fingerprints/strict DB, V3 ChemBERTa, V4 ChemBERTa-MTR).

Design (see project discussion):
  * Molecule pool  = FDA-approved drugs present in EVERY version's grouped database
                     (intersection pool) -> valid & retrievable under loose AND strict DBs.
  * Ground truth   = canonical "complete/loose" definition, decoupled from any version's
                     database: a molecule's targets from Chembl_36_filter.csv
                     (<=10000 nM, IC50/Ki/EC50, no relation/assay/lower-bound restriction).
  * Coverage parity= RDKit-parseable AND SMILES short enough that ChemBERTa does not
                     truncate (proxy: length <= MAX_SMILES_CHARS), so all four versions
                     operate on identical molecules/representations.
  * Frozen         = sampled once with a fixed seed and written to disk; every version
                     reads this exact file for both prediction and validation.

Output: query_master_FDA.csv  with columns [Ligand-id, SMILES, Target-id]
        (same 3-col format the existing pipelines already consume).
"""
from __future__ import annotations
import pandas as pd
from rdkit import Chem

CODE = "/Users/macbook/chembl/Code"
GH   = f"{CODE}/Github_target_prediction"
V2   = f"{CODE}/v2"

FILTER_CSV = f"{GH}/Chembl_36_filter.csv"          # canonical loose ground-truth source
V1_DB      = f"{GH}/Chembl_36_database.csv"         # grouped loose DB (Ligand-id, SMILES, Target-id)
V2_DB      = f"{V2}/v2_chembl_36_database.csv"      # grouped strict DB
FDA_CSV    = f"{GH}/FDA_data_initial.csv"           # all FDA-approved molecules
OUT        = f"{CODE}/query_master_FDA.csv"

N_QUERIES        = 100
RANDOM_STATE     = 42
MAX_SMILES_CHARS = 350   # conservative proxy for ChemBERTa max_length=384 tokens


def strip_chembl(s: pd.Series) -> pd.Series:
    return s.astype(str).str.replace("CHEMBL", "", regex=False)


def main() -> None:
    # 1) Canonical loose ground truth: molecule -> set(targets) from the filter file.
    print("[1] building canonical loose ground truth from Chembl_36_filter.csv ...")
    f = pd.read_csv(FILTER_CSV, usecols=["molecule_chembl_id", "target_chembl_id"], dtype=str)
    f["lig"] = strip_chembl(f["molecule_chembl_id"])
    f["tgt"] = strip_chembl(f["target_chembl_id"])
    gt = (f.groupby("lig")["tgt"]
            .agg(lambda s: ":".join(sorted(set(s), key=int)))
            .to_dict())
    print(f"    ground truth available for {len(gt):,} molecules")

    # 2) Canonical SMILES per ligand (from the loose grouped DB).
    v1 = pd.read_csv(V1_DB, dtype=str)
    v1["Ligand-id"] = v1["Ligand-id"].astype(str)
    smiles = dict(zip(v1["Ligand-id"], v1["SMILES"]))
    v1_ligs = set(v1["Ligand-id"])

    # 3) Intersection-eligible FDA pool: approved AND in both grouped databases.
    v2 = pd.read_csv(V2_DB, dtype=str)
    v2_ligs = set(v2["Ligand-id"].astype(str))
    fda = set(strip_chembl(pd.read_csv(FDA_CSV, dtype=str)["molecule_chembl_id"]))
    pool = fda & v1_ligs & v2_ligs
    print(f"[2] intersection-eligible FDA pool: {len(pool):,} "
          f"(FDA={len(fda):,}, v1db={len(v1_ligs):,}, v2db={len(v2_ligs):,})")

    # 4) Coverage parity: has ground truth + SMILES, RDKit-parseable, not ChemBERTa-truncated.
    eligible, drop_gt, drop_smi, drop_parse, drop_len = [], 0, 0, 0, 0
    for lig in pool:
        if lig not in gt:
            drop_gt += 1; continue
        smi = smiles.get(lig)
        if not isinstance(smi, str) or not smi:
            drop_smi += 1; continue
        if Chem.MolFromSmiles(smi) is None:
            drop_parse += 1; continue
        if len(smi) > MAX_SMILES_CHARS:
            drop_len += 1; continue
        eligible.append(lig)
    print(f"[3] coverage-parity filtered pool: {len(eligible):,} eligible "
          f"(dropped: no-GT={drop_gt}, no-SMILES={drop_smi}, unparseable={drop_parse}, "
          f"too-long={drop_len})")
    if len(eligible) < N_QUERIES:
        raise SystemExit(f"Only {len(eligible)} eligible molecules; need {N_QUERIES}.")

    # 5) Freeze: sample once with a fixed seed from a deterministically ordered pool.
    chosen = (pd.Series(sorted(eligible, key=int))
                .sample(n=N_QUERIES, random_state=RANDOM_STATE)
                .tolist())

    # 6) Write master file (3-col format the pipelines consume).
    out = pd.DataFrame({"Ligand-id": chosen,
                        "SMILES": [smiles[l] for l in chosen],
                        "Target-id": [gt[l] for l in chosen]})
    out.to_csv(OUT, index=False)
    npt = out["Target-id"].str.split(":").map(len)
    print(f"[4] wrote {OUT}  ({len(out)} queries)")
    print(f"    canonical true targets/query: mean={npt.mean():.1f} "
          f"min={npt.min()} max={npt.max()}")


if __name__ == "__main__":
    main()
