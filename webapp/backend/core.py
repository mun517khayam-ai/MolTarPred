# core.py — MolTarPred v2 prediction logic (framework-agnostic).
# Ported from the original Streamlit moltarpred.py; no Streamlit dependency.
# Loads precomputed artifacts produced by build_index.py.

from __future__ import annotations
import os
import pickle
from functools import lru_cache

import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
FPS_PATH = os.path.join(DATA_DIR, "fingerprints.pkl")
DB_PATH = os.path.join(DATA_DIR, "database.pkl")
TARGETS_PATH = os.path.join(DATA_DIR, "targets.pkl")
STATS_PATH = os.path.join(DATA_DIR, "stats.pkl")

CHEMBL_TARGET_URL = "https://www.ebi.ac.uk/chembl/explore/target/"
CHEMBL_COMPOUND_URL = "https://www.ebi.ac.uk/chembl/explore/compound/"

_morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


class Index:
    """Holds the in-memory database, fingerprints and target map. Built once."""

    def __init__(self):
        with open(FPS_PATH, "rb") as fh:
            self.fingerprints = pickle.load(fh)
        self.db = pd.read_pickle(DB_PATH)                 # ligand_id, smiles, target_ids
        targets = pd.read_pickle(TARGETS_PATH)            # target_chembl_id -> name, organism
        self.target_lookup = targets.set_index("target_chembl_id")
        with open(STATS_PATH, "rb") as fh:
            self.stats = pickle.load(fh)                  # written by build_index.py
        # Fast ligand_id -> target_ids map for supporting-hit lookups
        self._ligand_targets = dict(
            zip(self.db["ligand_id"].astype(str), self.db["target_ids"].astype(str))
        )

    @property
    def n_molecules(self) -> int:
        return len(self.db)

    @property
    def n_targets(self) -> int:
        return len(self.target_lookup)

    @property
    def n_interactions(self) -> int:
        return self.stats["n_interactions"]


@lru_cache(maxsize=1)
def get_index() -> Index:
    return Index()


def smiles_to_fingerprint(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return _morgan_gen.GetFingerprint(mol)


# Placeholder target names to drop *only when they have no organism*. A target
# with a real organism is kept even if its name is odd/unnamed.
_PLACEHOLDER_NAMES = {"unchecked"}


def _clean_org(value) -> str:
    s = str(value).strip()
    return "" if s.lower() in ("nan", "none", "") else s


def _target_meta(tid: str):
    chembl_id = f"CHEMBL{tid}"
    idx = get_index().target_lookup
    if chembl_id in idx.index:
        row = idx.loc[chembl_id]
        return chembl_id, str(row["target_name"]).strip(), _clean_org(row["target_organism"])
    return chembl_id, "Unknown", ""


def _is_placeholder(name: str, organism: str) -> bool:
    """A ChEMBL placeholder target (e.g. 'Unchecked') with no specific organism."""
    return name.lower() in _PLACEHOLDER_NAMES and not organism


def predict(query_smiles: str, top_k: int = 10) -> dict | None:
    """Run the full MolTarPred v2 prediction for one SMILES.

    Returns a dict with `targets` and `top_hits`, or None if the SMILES is invalid.
    """
    query_fp = smiles_to_fingerprint(query_smiles)
    if query_fp is None:
        return None

    ix = get_index()
    sims = DataStructs.BulkTanimotoSimilarity(query_fp, ix.fingerprints)

    db = ix.db.copy()
    db["similarity"] = sims
    top_hits = db.sort_values("similarity", ascending=False).head(top_k).reset_index(drop=True)

    # Count how many of the top-k hits support each target, and track which ones
    target_hit_counts: dict[str, int] = {}
    target_ligands: dict[str, list[str]] = {}
    for hit in top_hits.itertuples():
        ligand_chembl = f"CHEMBL{hit.ligand_id}"
        seen = set()
        for tid in str(hit.target_ids).split(":"):
            tid = tid.strip()
            if tid and tid not in seen:
                target_hit_counts[tid] = target_hit_counts.get(tid, 0) + 1
                target_ligands.setdefault(tid, []).append(ligand_chembl)
                seen.add(tid)

    targets = []
    for tid, count in target_hit_counts.items():
        chembl_id, name, organism = _target_meta(tid)
        if _is_placeholder(name, organism):
            continue
        targets.append({
            "target_chembl_id": chembl_id,
            "target_name": name,
            "organism": organism or "—",
            "reliability": round(count / top_k, 2),
            "supporting_hits": count,
            "supporting_ligands": target_ligands[tid],
            "chembl_url": CHEMBL_TARGET_URL + chembl_id,
        })
    targets.sort(key=lambda t: (t["reliability"], t["supporting_hits"]), reverse=True)

    hits = []
    for hit in top_hits.itertuples():
        chembl_id = f"CHEMBL{hit.ligand_id}"
        hits.append({
            "ligand_chembl_id": chembl_id,
            "smiles": hit.smiles,
            "similarity": round(float(hit.similarity), 4),
            "chembl_url": CHEMBL_COMPOUND_URL + chembl_id,
        })

    return {"targets": targets, "top_hits": hits}


def hits_for_target(target_chembl_id: str, top_hits: list[dict]) -> list[dict]:
    """From a list of top-hit dicts, return those active against the given target."""
    tid = target_chembl_id.replace("CHEMBL", "")
    ix = get_index()
    supporting = []
    for hit in top_hits:
        ligand_id = hit["ligand_chembl_id"].replace("CHEMBL", "")
        target_ids = ix._ligand_targets.get(ligand_id, "").split(":")
        if tid in target_ids:
            supporting.append(hit)
    return supporting
