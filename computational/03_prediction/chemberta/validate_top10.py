#!/usr/bin/env python3
"""
Validation, faithfully replicating Top10.py's confusion-matrix logic, plus F1.

For each query, compares predicted targets (from a *-Top10 csv) against the true
targets (query_10nM.csv). TN universe = all unique targets seen across true+pred
(exactly as Top10.py defines `total`). Reports per-query metrics and the average.

Usage:
  python validate_top10.py --pred ../ChemBERTa-Top10_raw.csv --method ChemBERTa_raw
  python validate_top10.py --pred /.../Github_target_prediction/MolTarPred-Top10.csv --method MolTarPred
"""
from __future__ import annotations
import argparse, csv, math
from pathlib import Path
import numpy as np

GH = "/Users/macbook/chembl/Code/Github_target_prediction"
# v2 strict setup: validate against the master query set with strict ground truth
V2_QUERY = "/Users/macbook/chembl/Code/query_master_FDA_strictGT.csv"


def read_top(path):
    rows = []
    with open(path) as f:
        for r in csv.reader(f):
            rows.append((r[0], r[1], r[2].split(":")))
    return rows[1:]  # drop header


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred", required=True)
    ap.add_argument("--query", default=V2_QUERY)
    ap.add_argument("--method", default="method")
    ap.add_argument("--outdir", default=str(Path(__file__).parent))
    args = ap.parse_args()

    db1 = read_top(args.pred)          # predictions: (id, smiles, [pred targets])
    qr1 = read_top(args.query)         # truth:       (id, smiles, [true targets])

    # build TN universe + aligned (name, true, pred) triples — exactly as Top10.py
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

    rows = []  # per-query
    for name, real, pred in val_rel:
        TP = sum(1 for k in pred if k in real)
        FP = len(pred) - TP
        FN = len(real) - TP
        TN = total - FN - TP - FP
        acc = (TP + TN) / (TP + FP + TN + FN)
        ppv = TP / (TP + FP) if (TP + FP) else 0.0
        recall = TP / (TP + FN) if (TP + FN) else 0.0
        f1 = 2 * ppv * recall / (ppv + recall) if (ppv + recall) else 0.0
        den = math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
        mcc = (TP*TN - FP*FN) / den if den else 0.0
        rows.append(dict(lig=name, NPT=len(pred), acc=acc, ppv=ppv, recall=recall,
                         f1=f1, mcc=mcc, TN=TN, FP=FP, FN=FN, TP=TP))

    def avg(k): return float(np.mean([r[k] for r in rows]))
    summary = dict(Method=args.method, avNPT=avg("NPT"), avAccuracy=avg("acc"),
                   avPrecision=avg("ppv"), avRecall=avg("recall"), avF1=avg("f1"),
                   avMCC=avg("mcc"), avTN=avg("TN"), avFP=avg("FP"), avFN=avg("FN"), avTP=avg("TP"))

    outdir = Path(args.outdir)
    per = outdir / f"{args.method}_per.csv"
    with open(per, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Ligand-ID", "NPT", "Accuracy", "Precision", "Recall", "F1", "MCC", "TN", "FP", "FN", "TP"])
        for r in rows:
            w.writerow([r["lig"], r["NPT"], r["acc"], r["ppv"], r["recall"], r["f1"], r["mcc"], r["TN"], r["FP"], r["FN"], r["TP"]])
    ave = outdir / f"Validation_{args.method}_ave.csv"
    with open(ave, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(list(summary.keys()))
        w.writerow([f"{v:.4f}" if isinstance(v, float) else v for v in summary.values()])

    print(f"[{args.method}]  NPT {summary['avNPT']:.1f}  Recall {summary['avRecall']:.3f}  "
          f"Precision {summary['avPrecision']:.3f}  F1 {summary['avF1']:.3f}  MCC {summary['avMCC']:.3f}")
    print(f"  wrote {ave.name} and {per.name}")


if __name__ == "__main__":
    main()
