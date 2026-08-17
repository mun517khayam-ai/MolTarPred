#!/usr/bin/env python3
"""
V1 (fingerprints / LOOSE database) run on the shared master query set.

Self-contained: does not touch any existing V1 script or file. Reads the frozen
master query set + the V1 loose grouped database, masks the master queries out of
the reference library, runs canonical MolTarPred (Morgan r2/2048 + Tanimoto, top-10,
union of neighbour targets), and validates against the master canonical ground truth.

All outputs are written into this master_run/ subfolder.
"""
from __future__ import annotations
import csv, math
from pathlib import Path
from time import time
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

# ---- version-specific config -------------------------------------------------
CODE     = "/Users/macbook/chembl/Code"
DB_CSV   = f"{CODE}/Github_target_prediction/Chembl_36_database.csv"   # LOOSE grouped DB
QUERY    = f"{CODE}/query_master_FDA_strictGT.csv"                              # shared master queries
OUTDIR   = Path(__file__).parent
METHOD   = "V1_Morgan_Tanimoto_strictGT"
TOPK     = 10
# -----------------------------------------------------------------------------


def load_predictions():
    t0 = time()
    q = pd.read_csv(QUERY, dtype=str)
    qset = set(q["Ligand-id"].astype(str))

    db = pd.read_csv(DB_CSV, dtype=str)
    db["Ligand-id"] = db["Ligand-id"].astype(str)
    before = len(db)
    db = db[~db["Ligand-id"].isin(qset)]                    # mask master queries out
    print(f"[db] {DB_CSV.split('/')[-1]}: masked {before-len(db)} query ligands "
          f"({before} -> {len(db)})", flush=True)

    db["Mol"] = db["SMILES"].apply(Chem.MolFromSmiles)
    db = db.dropna(subset=["Mol"]).reset_index(drop=True)
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    db_fps = [gen.GetFingerprint(m) for m in db["Mol"]]
    db_tgt = db["Target-id"].astype(str).tolist()
    print(f"[db] {len(db_fps)} fingerprints built ({time()-t0:.0f}s)", flush=True)

    q["Mol"] = q["SMILES"].apply(Chem.MolFromSmiles)
    q = q.dropna(subset=["Mol"]).reset_index(drop=True)

    rows = []
    for _, qr in q.iterrows():
        qfp = gen.GetFingerprint(qr["Mol"])
        sims = DataStructs.BulkTanimotoSimilarity(qfp, db_fps)
        order = sorted(range(len(sims)), key=lambda i: sims[i], reverse=True)[:TOPK]
        uniq = set()
        for i in order:
            uniq.update(db_tgt[i].split(":"))
        uniq.discard("")
        rows.append((qr["Ligand-id"], qr["SMILES"], ":".join(sorted(uniq, key=int))))

    pred_path = OUTDIR / f"{METHOD}_Top{TOPK}.csv"
    pd.DataFrame(rows, columns=["Ligand-ID", "SMILES", "Targets-ID"]).to_csv(pred_path, index=False)
    npt = np.mean([len(r[2].split(":")) if r[2] else 0 for r in rows])
    print(f"[predict] wrote {pred_path.name}  avg predicted targets/query={npt:.1f} "
          f"({time()-t0:.0f}s)", flush=True)
    return pred_path


def validate(pred_path):
    def read_top(path):
        with open(path) as f:
            return [(r[0], r[1], r[2].split(":")) for r in csv.reader(f)][1:]

    db1 = read_top(pred_path)          # predictions
    qr1 = read_top(QUERY)              # truth (master ground truth)

    uniq_targets, val_rel = [], []
    for lig in qr1:
        ytnad, qrtarg = lig[0], lig[2]
        for x in db1:
            if x[0] == ytnad:
                outarg = x[2]
                val_rel.append((ytnad, qrtarg, outarg))
                for t in outarg + qrtarg:
                    if t not in uniq_targets:
                        uniq_targets.append(t)
    total = len(uniq_targets)

    rows = []
    for name, real, pred in val_rel:
        TP = sum(1 for k in pred if k in real)
        FP = len(pred) - TP
        FN = len(real) - TP
        TN = total - FN - TP - FP
        acc = (TP + TN) / (TP + FP + TN + FN) if (TP + FP + TN + FN) else 0.0
        ppv = TP / (TP + FP) if (TP + FP) else 0.0
        rec = TP / (TP + FN) if (TP + FN) else 0.0
        f1 = 2 * ppv * rec / (ppv + rec) if (ppv + rec) else 0.0
        den = math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
        mcc = (TP*TN - FP*FN) / den if den else 0.0
        rows.append(dict(lig=name, NPT=len(pred), acc=acc, ppv=ppv, recall=rec,
                         f1=f1, mcc=mcc, TN=TN, FP=FP, FN=FN, TP=TP))

    def avg(k): return float(np.mean([r[k] for r in rows]))
    summary = dict(Method=METHOD, avNPT=avg("NPT"), avAccuracy=avg("acc"),
                   avPrecision=avg("ppv"), avRecall=avg("recall"), avF1=avg("f1"),
                   avMCC=avg("mcc"), avTN=avg("TN"), avFP=avg("FP"), avFN=avg("FN"), avTP=avg("TP"))

    with open(OUTDIR / f"{METHOD}_per.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Ligand-ID","NPT","Accuracy","Precision","Recall","F1","MCC","TN","FP","FN","TP"])
        for r in rows:
            w.writerow([r["lig"],r["NPT"],r["acc"],r["ppv"],r["recall"],r["f1"],r["mcc"],r["TN"],r["FP"],r["FN"],r["TP"]])
    with open(OUTDIR / f"Validation_{METHOD}_ave.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(list(summary.keys()))
        w.writerow([f"{v:.4f}" if isinstance(v, float) else v for v in summary.values()])

    print(f"[{METHOD}]  n={len(rows)}  NPT {summary['avNPT']:.1f}  "
          f"Recall {summary['avRecall']:.3f}  Precision {summary['avPrecision']:.3f}  "
          f"F1 {summary['avF1']:.3f}  MCC {summary['avMCC']:.3f}", flush=True)


if __name__ == "__main__":
    validate(load_predictions())
