# build_index.py — one-off precompute of the MolTarPred v2 search index.
# Produces backend/data/{fingerprints.pkl, database.pkl, targets.pkl, stats.pkl} so
# the API starts instantly instead of fingerprinting ~900k molecules on every boot.
#
# Re-run this whenever the v2 database is rebuilt; see Code/DATA_PIPELINE_FIXES.md.
#
# Usage:  python build_index.py

import os
import pickle
import time

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")

# --- v2 sources (MolTarPred v2) ---------------------------------------------
# NB: this is the UNMASKED database, and that is deliberate.
# The *_masked.csv variant has the 100 benchmark query ligands removed, which is
# required for validation (a query must not retrieve itself) but wrong for the
# deployed tool: those are approved drugs, and withholding them would silently
# impoverish the neighbourhood returned for any user searching near that chemistry.
# Masking belongs to the evaluation pipeline only -- do not point this at
# *_masked.csv to "match" the reported metrics.
V2_DB = "/Users/macbook/chembl/Code/v2/v2_chembl_36_database.csv"
V2_FILTER = "/Users/macbook/chembl/Code/v2/v2_chembl_36_filter.csv"

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
os.makedirs(DATA_DIR, exist_ok=True)

gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def main():
    t0 = time.time()

    db = pd.read_csv(V2_DB, dtype=str)
    db.columns = ["ligand_id", "smiles", "target_ids"]
    print(f"Loaded v2 database: {len(db):,} rows")

    fps, valid = [], []
    for idx, smi in enumerate(db["smiles"]):
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            fps.append(gen.GetFingerprint(mol))
            valid.append(idx)
    db_valid = db.iloc[valid].reset_index(drop=True)
    print(f"Fingerprinted {len(fps):,} valid molecules "
          f"({len(db) - len(fps)} failed)  [{time.time()-t0:.0f}s]")

    # Target ids actually reachable by a search, i.e. present on a fingerprinted
    # molecule. The filter file also contains targets carried only by molecules that
    # never reach the index (unparseable SMILES), so counting it directly would
    # overstate what is predictable.
    db_target_ids = {
        t for row in db_valid["target_ids"].astype(str) for t in row.split(":") if t
    }
    n_interactions = sum(
        len({t for t in row.split(":") if t}) for row in db_valid["target_ids"].astype(str)
    )

    # Slim target map: one row per target with name + organism
    filt = pd.read_csv(V2_FILTER, usecols=["target_chembl_id", "target_name", "target_organism"])
    targets = filt.drop_duplicates("target_chembl_id").reset_index(drop=True)
    targets = targets[
        targets["target_chembl_id"].str.replace("CHEMBL", "", regex=False).isin(db_target_ids)
    ].reset_index(drop=True)
    print(f"Target map: {len(targets):,} targets present in the searchable database "
          f"({len(db_target_ids) - len(targets)} db targets missing from the filter map)")
    print(f"Ligand-target interactions: {n_interactions:,}")

    stats = {
        "n_molecules": len(db_valid),
        "n_targets": len(targets),
        "n_interactions": n_interactions,
    }

    with open(os.path.join(DATA_DIR, "fingerprints.pkl"), "wb") as fh:
        pickle.dump(fps, fh, protocol=pickle.HIGHEST_PROTOCOL)
    db_valid.to_pickle(os.path.join(DATA_DIR, "database.pkl"))
    targets.to_pickle(os.path.join(DATA_DIR, "targets.pkl"))
    with open(os.path.join(DATA_DIR, "stats.pkl"), "wb") as fh:
        pickle.dump(stats, fh, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"Artifacts written to {DATA_DIR}  [{time.time()-t0:.0f}s total]")


if __name__ == "__main__":
    main()
