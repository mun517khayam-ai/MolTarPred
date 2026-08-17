#!/usr/bin/env python3
"""
Fusion of ChemBERTa (whiten-64 cosine) and Morgan/Tanimoto for MolTarPred-style
target prediction. Same DB + queries + top-10 union logic; only the *ranking* is fused.

Stage 1 (cached): per-query similarity matrices for both metrics over all DB molecules
  artifacts/cb_sim.npy   [Q, Ndb] ChemBERTa whiten-64 cosine
  artifacts/tan_sim.npy  [Q, Ndb] Morgan/Tanimoto
Stage 2: fuse and write a *-Top10 csv (RRF or weighted-score).

Usage:
  python fusion_MolTarPred.py --method rrf --rrf-k 60
  python fusion_MolTarPred.py --method score --alpha 0.5
"""
from __future__ import annotations
import argparse, csv
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

ART = Path(__file__).parent / "artifacts"
GH = "/Users/macbook/chembl/Code/Github_target_prediction"


def fit_whitening(X, dim, eps=1e-6):
    mu = X.mean(0); Xc = X - mu
    vals, vecs = np.linalg.eigh((Xc.T @ Xc) / (len(X) - 1))
    o = np.argsort(vals)[::-1]; vals, vecs = vals[o][:dim], vecs[:, o][:, :dim]
    return mu, vecs / np.sqrt(vals + eps)


def l2(X):
    return X / np.linalg.norm(X, axis=1, keepdims=True).clip(min=1e-12)


def build_matrices(art, whiten_dim=64):
    """Compute + cache the two [Q, Ndb] similarity matrices.

    NB: the cache is reused WITHOUT validation. cb_sim/tan_sim are keyed to the
    database they were built from, so after any database rebuild they must be
    DELETED, not merely regenerated, or fusion reads misaligned rows.
    """
    cb_path, tan_path = art/"cb_sim.npy", art/"tan_sim.npy"
    if cb_path.exists() and tan_path.exists():
        print("[cache] loading cb_sim + tan_sim")
        return np.load(cb_path, mmap_mode="r"), np.load(tan_path, mmap_mode="r")

    db = pd.read_parquet(art/"db_meta.parquet"); qr = pd.read_parquet(art/"query_meta.parquet")
    demb = np.load(art/"db_emb.npy").astype(np.float64)
    qemb = np.load(art/"query_emb.npy").astype(np.float64)
    mu, W = fit_whitening(demb, whiten_dim)
    demb_w, qemb_w = l2((demb-mu)@W), l2((qemb-mu)@W)
    print("[build] ChemBERTa cosine matrix")
    cb = (qemb_w @ demb_w.T).astype(np.float32)
    np.save(cb_path, cb)

    print("[build] Morgan fingerprints (DB + query)")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    db_fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in db["smiles"]]
    q_fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in qr["smiles"]]
    print("[build] Tanimoto matrix")
    tan = np.empty((len(q_fps), len(db_fps)), dtype=np.float32)
    for i, qf in enumerate(q_fps):
        tan[i] = DataStructs.BulkTanimotoSimilarity(qf, db_fps)
    np.save(tan_path, tan)
    return cb, tan


def ranks_desc(sim_row):
    """rank 0 = most similar (largest sim)."""
    order = np.argsort(-sim_row, kind="stable")
    r = np.empty_like(order); r[order] = np.arange(len(order))
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", choices=["rrf", "score"], default="rrf")
    ap.add_argument("--rrf-k", type=int, default=60)
    ap.add_argument("--alpha", type=float, default=0.5, help="score fusion: alpha*cb + (1-alpha)*tanimoto")
    ap.add_argument("--topk", type=int, default=10)
    ap.add_argument("--artdir", default=str(ART),
                    help="embedding artifacts dir; one per backbone (see run_chemberta.sh)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    art = Path(args.artdir)
    cb, tan = build_matrices(art)
    db = pd.read_parquet(art/"db_meta.parquet"); qr = pd.read_parquet(art/"query_meta.parquet")
    db_targets = [t.split(":") for t in db["targets"].astype(str)]

    tag = f"rrf{args.rrf_k}" if args.method == "rrf" else f"score{args.alpha}"
    rows = []
    for qi in range(len(qr)):
        cb_row, tan_row = np.asarray(cb[qi]), np.asarray(tan[qi])
        if args.method == "rrf":
            fused = 1.0/(args.rrf_k + ranks_desc(cb_row)) + 1.0/(args.rrf_k + ranks_desc(tan_row))
        else:  # min-max normalise each then weight
            def mm(v): return (v - v.min())/(v.max()-v.min()+1e-12)
            fused = args.alpha*mm(cb_row) + (1-args.alpha)*mm(tan_row)
        top = np.argpartition(-fused, args.topk)[:args.topk]
        top = top[np.argsort(-fused[top])]
        uniq = set()
        for d in top:
            uniq.update(db_targets[d])
        uniq.discard("")
        rows.append((qr["ligand_id"].iloc[qi], qr["smiles"].iloc[qi], ":".join(sorted(uniq, key=int))))

    out = Path(args.out) if args.out else art.parent / f"Fusion-Top{args.topk}_{tag}.csv"
    pd.DataFrame(rows, columns=["Ligand-ID", "SMILES", "Targets-ID"]).to_csv(out, index=False)
    npt = np.mean([len(r[2].split(":")) if r[2] else 0 for r in rows])
    print(f"[done] wrote {out.name}  | method={tag}  avg predicted targets/query={npt:.1f}")


if __name__ == "__main__":
    main()
