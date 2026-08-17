#!/usr/bin/env python3
"""
ChemBERTa analogue of Initial_MolTarPred.py.

Same workflow as MolTarPred (top-10 nearest neighbours -> union of their targets),
but similarity = ChemBERTa cosine instead of Morgan/Tanimoto. Reads the cached
embeddings from embed_db.py.

Output (identical format to MolTarPred-Top10.csv):
  Ligand-ID, SMILES, Targets-ID   (Targets-ID = ':'-joined, sorted ascending by int)
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np, pandas as pd

ART = Path(__file__).parent / "artifacts"


def fit_whitening(X, dim, eps=1e-6):
    mu = X.mean(0)
    Xc = X - mu
    vals, vecs = np.linalg.eigh((Xc.T @ Xc) / (len(X) - 1))
    o = np.argsort(vals)[::-1]
    vals, vecs = vals[o][:dim], vecs[:, o][:, :dim]
    return mu, vecs / np.sqrt(vals + eps)


def l2(X):
    return X / np.linalg.norm(X, axis=1, keepdims=True).clip(min=1e-12)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--artdir", default=str(ART))
    ap.add_argument("--topk", type=int, default=10)
    ap.add_argument("--whiten-dim", type=int, default=0, help="0 = raw cosine; else whiten to N dims")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    art = Path(args.artdir)

    db = pd.read_parquet(art / "db_meta.parquet")
    qr = pd.read_parquet(art / "query_meta.parquet")
    demb = np.load(art / "db_emb.npy").astype(np.float64)
    qemb = np.load(art / "query_emb.npy").astype(np.float64)

    tag = "raw"
    if args.whiten_dim:
        mu, W = fit_whitening(demb, args.whiten_dim)      # fit on the DB
        demb = l2((demb - mu) @ W)
        qemb = l2((qemb - mu) @ W)                        # SAME transform for queries
        tag = f"whiten{args.whiten_dim}"
    else:
        demb, qemb = l2(demb), l2(qemb)                   # already normalized; idempotent

    db_targets = [t.split(":") for t in db["targets"].astype(str)]
    sims = qemb @ demb.T                                  # [Q, Ndb] cosine
    out_rows = []
    for qi in range(len(qr)):
        top = np.argpartition(-sims[qi], args.topk)[:args.topk]
        top = top[np.argsort(-sims[qi][top])]             # ordered top-k hits
        uniq = set()
        for d in top:
            uniq.update(db_targets[d])
        uniq.discard("")
        targ = ":".join(sorted(uniq, key=int))
        out_rows.append((qr["ligand_id"].iloc[qi], qr["smiles"].iloc[qi], targ))

    out_path = Path(args.out) if args.out else art.parent / f"ChemBERTa-Top{args.topk}_{tag}.csv"
    pd.DataFrame(out_rows, columns=["Ligand-ID", "SMILES", "Targets-ID"]).to_csv(out_path, index=False)
    npt = np.mean([len(r[2].split(":")) if r[2] else 0 for r in out_rows])
    print(f"[done] wrote {out_path}  | space={tag}  avg predicted targets/query={npt:.1f}")


if __name__ == "__main__":
    main()
