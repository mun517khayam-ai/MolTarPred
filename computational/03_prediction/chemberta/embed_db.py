#!/usr/bin/env python3
"""
Embed the MolTarPred database + query molecules with ChemBERTa,
masked-mean-pool + L2-normalize (so a dot product is a cosine similarity).

Both backbones reported in the study run through this one script; they differ
only in --model:
  seyonec/ChemBERTa-zinc-base-v1   MLM-pretrained on ZINC, 768-dim, 44.1M params
  DeepChem/ChemBERTa-77M-MTR       multi-task regression over 199 properties,
                                   384-dim, 3.4M params
Give each its own --outdir; the artifacts are keyed to the model that wrote them.

Molecules whose SMILES fail RDKit parsing are dropped, so the candidate pool
matches the fingerprint pipeline's.

Outputs (in --outdir):
  db_emb.npy        float32 [Nvalid, hidden]   row i <-> db_meta row i
  db_meta.parquet   columns: [ligand_id, smiles, targets]   (targets = ':'-joined)
  query_emb.npy     float32 [Qvalid, hidden]
  query_meta.parquet
  meta.json
"""
from __future__ import annotations
import argparse, json, time
from pathlib import Path
import numpy as np, pandas as pd
import torch, torch.nn.functional as F
from transformers import AutoTokenizer, AutoModel
from rdkit import Chem

MODEL_NAME = "seyonec/ChemBERTa-zinc-base-v1"   # default; override with --model
CODE = "/Users/macbook/chembl/Code"
# v2 strict setup: stricter-filtered database (master queries removed) + master query set (strict ground truth)
V2_DB = f"{CODE}/v2/v2_chembl_36_database_masked.csv"
V2_QUERY = f"{CODE}/query_master_FDA_strictGT.csv"


def device_auto():
    return "mps" if torch.backends.mps.is_available() else ("cuda" if torch.cuda.is_available() else "cpu")


def load_valid(csv_path):
    """Read a MolTarPred 3-col csv; keep rows with RDKit-parseable SMILES."""
    df = pd.read_csv(csv_path, dtype=str)
    df.columns = ["ligand_id", "smiles", "targets"]
    keep = [Chem.MolFromSmiles(s) is not None for s in df["smiles"]]
    n_drop = len(df) - int(np.sum(keep))
    df = df[keep].reset_index(drop=True)
    print(f"[load] {csv_path.split('/')[-1]}: {len(df)} valid ({n_drop} dropped)", flush=True)
    return df


@torch.no_grad()
def embed_frame(df, tok, mdl, device, bs, max_len):
    smiles = df["smiles"].to_numpy()
    N = len(smiles)
    out = np.empty((N, mdl.config.hidden_size), dtype=np.float32)
    order = np.argsort(-np.char.str_len(smiles.astype(str)), kind="stable")  # long->short
    t0 = time.time()
    for s in range(0, N, bs):
        sel = order[s:s+bs]
        enc = tok(list(smiles[sel]), padding=True, truncation=True,
                  max_length=max_len, return_tensors="pt").to(device)
        h = mdl(**enc).last_hidden_state
        m = enc["attention_mask"].unsqueeze(-1).to(h.dtype)
        pooled = (h*m).sum(1) / m.sum(1).clamp(min=1e-9)
        out[sel] = F.normalize(pooled, dim=-1).float().cpu().numpy()
        if s % (bs*50) == 0 or s+bs >= N:
            done = min(s+bs, N); rate = done/max(time.time()-t0, 1e-9)
            print(f"    {done}/{N} ({100*done/N:5.1f}%) {rate:6.1f} mol/s ETA {(N-done)/max(rate,1e-9)/60:5.1f}m", flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=MODEL_NAME,
                    help="HuggingFace checkpoint; see the module docstring for the two used here")
    ap.add_argument("--db", default=V2_DB)
    ap.add_argument("--query", default=V2_QUERY)
    ap.add_argument("--outdir", default=str(Path(__file__).parent / "artifacts"))
    ap.add_argument("--batch-size", type=int, default=96)
    ap.add_argument("--max-length", type=int, default=384)
    args = ap.parse_args()
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    dev = device_auto()

    tok = AutoTokenizer.from_pretrained(args.model)
    mdl = AutoModel.from_pretrained(args.model).to(dev).eval()
    print(f"[model] {args.model} on {dev}", flush=True)

    qdf = load_valid(args.query)
    print("[embed] query molecules", flush=True)
    qemb = embed_frame(qdf, tok, mdl, dev, args.batch_size, args.max_length)
    np.save(out/"query_emb.npy", qemb); qdf.to_parquet(out/"query_meta.parquet", index=False)

    ddf = load_valid(args.db)
    print("[embed] database molecules", flush=True)
    demb = embed_frame(ddf, tok, mdl, dev, args.batch_size, args.max_length)
    np.save(out/"db_emb.npy", demb); ddf.to_parquet(out/"db_meta.parquet", index=False)

    (out/"meta.json").write_text(json.dumps({
        "model": args.model, "pooling": "masked_mean+l2norm",
        "n_db": len(ddf), "n_query": len(qdf), "hidden": int(mdl.config.hidden_size),
        "max_length": args.max_length, "device": dev}, indent=2))
    print(f"[done] db {demb.shape} | query {qemb.shape}", flush=True)


if __name__ == "__main__":
    main()
