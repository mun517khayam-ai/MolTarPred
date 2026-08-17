#!/usr/bin/env bash
# Run the full ChemBERTa target-prediction pipeline for one backbone against the
# v2.2 strict database and the frozen strict-ground-truth benchmark.
#
#   ./run_chemberta.sh base   # seyonec/ChemBERTa-zinc-base-v1
#   ./run_chemberta.sh mtr    # DeepChem/ChemBERTa-77M-MTR
#
# Each backbone writes into its own artifacts/ and results/ subdirectory, since
# the cached embeddings and similarity matrices are keyed to the model.
set -euo pipefail
export KMP_DUPLICATE_LIB_OK=TRUE
cd "$(dirname "$0")"

case "${1:-}" in
  base) MODEL="seyonec/ChemBERTa-zinc-base-v1"; TAG="base" ;;
  mtr)  MODEL="DeepChem/ChemBERTa-77M-MTR";     TAG="mtr"  ;;
  *)    echo "usage: $0 {base|mtr}" >&2; exit 2 ;;
esac

ART="$(pwd)/artifacts_${TAG}"
OUT="$(pwd)/results_${TAG}"
mkdir -p "$OUT"

# Absolute paths as used in development -- edit these two for another machine.
DB="/Users/macbook/chembl/Code/v2/v2_chembl_36_database_masked.csv"
Q="/Users/macbook/chembl/Code/query_master_FDA_strictGT.csv"

echo "== [1/4] embed database + queries with $MODEL =="
# NB: delete artifacts_${TAG}/{cb_sim,tan_sim}.npy after any database rebuild --
# they are keyed to the database and are reused WITHOUT validation.
python3 embed_db.py --model "$MODEL" --db "$DB" --query "$Q" --outdir "$ART"

echo "== [2/4] ChemBERTa predictions (raw + whitened to 64 dims) =="
python3 ChemBERTa_MolTarPred.py --artdir "$ART" --whiten-dim 0 \
        --out "$OUT/ChemBERTa-Top10_raw.csv"
python3 ChemBERTa_MolTarPred.py --artdir "$ART" --whiten-dim 64 \
        --out "$OUT/ChemBERTa-Top10_whiten64.csv"

echo "== [3/4] fusion with Morgan/Tanimoto =="
for k in 10 30 60; do
  python3 fusion_MolTarPred.py --method rrf --rrf-k "$k" --artdir "$ART" \
          --out "$OUT/Fusion-Top10_rrf${k}.csv"
done
for a in 0.0 0.1 0.2 0.3 0.5 0.7; do
  python3 fusion_MolTarPred.py --method score --alpha "$a" --artdir "$ART" \
          --out "$OUT/Fusion-Top10_score${a}.csv"
done

echo "== [4/4] validation against the strict ground truth =="
python3 validate_top10.py --pred "$OUT/ChemBERTa-Top10_raw.csv"      --query "$Q" \
        --method ChemBERTa_raw      --outdir "$OUT"
python3 validate_top10.py --pred "$OUT/ChemBERTa-Top10_whiten64.csv" --query "$Q" \
        --method ChemBERTa_whiten64 --outdir "$OUT"
for k in 10 30 60; do
  python3 validate_top10.py --pred "$OUT/Fusion-Top10_rrf${k}.csv" --query "$Q" \
          --method "Fusion_rrf${k}" --outdir "$OUT"
done
for a in 0.0 0.1 0.2 0.3 0.5 0.7; do
  python3 validate_top10.py --pred "$OUT/Fusion-Top10_score${a}.csv" --query "$Q" \
          --method "Fusion_score${a}" --outdir "$OUT"
done

echo "== DONE ($TAG) -- metrics in $OUT =="
