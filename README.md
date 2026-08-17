# MolTarPred

**Ligand-based protein target prediction on ChEMBL 36.**

Predicts which protein targets a small molecule is likely to bind, by finding its
nearest neighbours in a curated bioactivity database and pooling their known targets.
Each prediction carries a **reliability score** — the fraction of the retrieved
neighbours annotated with that target.

This repository contains the computational pipeline and the deployed web application
for the MSc research project *"MolTarPred: A webserver to predict molecular targets of
drug molecules"* (Imperial College London, 2026), which re-implements and extends
MolTarPred (He, Caba & Ballester, *Digital Discovery*, 2025) on ChEMBL 36.

> **Central finding.** For ligand-based target prediction the molecular *representation*
> is not the bottleneck — **database curation and honest evaluation are.** A Morgan
> fingerprint with Tanimoto similarity on a well-curated database is a strong baseline;
> off-the-shelf chemical language models do not beat it; fusing the two reaches parity.
> What governs whether a prediction can be trusted is how close the query sits to known
> chemistry.

---

## Results

Mean over the 100-drug benchmark, strict ground truth, top-*k* = 10. Full per-query
distributions, statistics and figures are in the dissertation.

| Method | NPT | Precision | Recall | F1 | MCC |
|---|---:|---:|---:|---:|---:|
| MolTarPred v2.1 (loose filters) | 17.0 | 0.366 | 0.511 | 0.332 | 0.369 |
| **MolTarPred v2.2 (strict filters)** | 11.1 | 0.445 | 0.561 | 0.385 | 0.427 |
| ChemBERTa-base, whitened | 9.8 | 0.413 | 0.483 | 0.345 | 0.381 |
| ChemBERTa-base + Morgan fusion | 11.3 | 0.451 | 0.558 | 0.385 | 0.427 |
| ChemBERTa-MTR, whitened | 10.7 | 0.433 | 0.533 | 0.376 | 0.414 |
| **ChemBERTa-MTR + Morgan fusion** | 11.1 | 0.457 | 0.578 | 0.396 | 0.439 |


---

## Repository layout

```
computational/
├── 01_database/          ChEMBL 36 -> filtered, grouped ligand-target databases
│   ├── v2.1_chembl_36_filter.py     base filters (IC50/Ki/EC50, <=10,000 nM, >=10 nM)
│   ├── v2.1_chembl_36_database.py   group to one row per ligand
│   ├── v2.2_chembl_36_filter.py     adds relation '=', binding assays only, >0 nM
│   └── v2.2_chembl_36_database.py
├── 02_benchmark/         the frozen, version-independent benchmark
│   ├── build_master_query.py           intersection pool -> 100 drugs (seed 42)
│   ├── build_master_query_strictGT.py  strict ground-truth answer key
│   ├── exclude_queries_from_db.py      masks query ligands out of a database
│   ├── query_master_FDA.csv            100 molecules + loose ground truth
│   └── query_master_FDA_strictGT.csv   100 molecules + strict ground truth  <- the benchmark
└── 03_prediction/
    ├── morgan/           Morgan/Tanimoto k-NN for v2.1 and v2.2 (self-contained;
    │                     each script predicts and validates in one pass)
    └── chemberta/        both transformer backbones share one pipeline
        ├── embed_db.py             SMILES -> masked-mean-pool -> L2-normalised vectors
        ├── ChemBERTa_MolTarPred.py cosine k-NN, raw or whitened
        ├── fusion_MolTarPred.py    RRF / weighted-score fusion with Morgan
        ├── validate_top10.py       per-query confusion matrix and metrics
        └── run_chemberta.sh        driver: `./run_chemberta.sh {base|mtr}`

webapp/                   FastAPI + React application serving MolTarPred v2.2
├── backend/              main.py (API), core.py (prediction), build_index.py
└── frontend/             Vite + React single-page app
```

The two ChemBERTa backbones were originally maintained as separate directories whose
files were byte-identical apart from the model name. They are consolidated here into one
pipeline selected by `--model`; each backbone writes to its own `artifacts_*/` and
`results_*/` so the cached embeddings and similarity matrices stay keyed to their model.

This repository holds only what is needed to run the pipeline. Figures, the
plotting/statistics notebook, and the per-query result files are not included — the
results are reported in the dissertation.

---

## ⚠️ Paths are absolute — read before running

**Every script contains absolute paths under `/Users/macbook/chembl/`.** They were written
for a single development machine and are preserved here unchanged so that the code
matches exactly what produced the reported results.

To run any of this elsewhere you must edit the path constants at the top of each script.
They are declared as module-level constants (`CODE`, `DB_CSV`, `QUERY`, `V2_DB`, `ART`, …)
in the first ~30 lines of each file, so the edit is mechanical. The database connection
string is likewise hardcoded as `postgresql://macbook@localhost:5432/chembl_36` in the two
filter scripts.

This is a known limitation rather than an oversight, and it is the main thing to fix before
anyone else runs the pipeline.

---

## Prerequisites

- **Python 3.11** with the packages in `requirements.txt`
- **PostgreSQL** with a ChEMBL 36 dump restored locally
  ([download](https://chembl.gitbook.io/chembl-interface-documentation/downloads))
- **Node 18+** for the web app frontend
- Apple-Silicon GPU (MPS), CUDA, or CPU for the embedding step

```bash
conda create -n moltarpred python=3.11 && conda activate moltarpred
pip install -r requirements.txt
```

---

## Running the pipeline

Steps 1–3 take a few minutes; step 5 is the expensive one (~30 min of embedding).

```bash
# 1. Build both filtered databases from the local ChEMBL 36 instance   (~25 s each)
python computational/01_database/v2.1_chembl_36_filter.py
python computational/01_database/v2.1_chembl_36_database.py
python computational/01_database/v2.2_chembl_36_filter.py
python computational/01_database/v2.2_chembl_36_database.py

# 2. Mask the 100 benchmark ligands out of each reference database
python computational/02_benchmark/exclude_queries_from_db.py --db <grouped_db.csv>

# 3. Derive the strict ground-truth answer key for the frozen query set
python computational/02_benchmark/build_master_query_strictGT.py

# 4. Morgan/Tanimoto predictions + validation
python computational/03_prediction/morgan/run_master_v2_strictGT.py
python computational/03_prediction/morgan/run_master_v1_strictGT.py

# 5. ChemBERTa embeddings, predictions, fusion and validation (one call per backbone)
#    DELETE THE CACHED SIMILARITY MATRICES FIRST after any database rebuild (see below)
cd computational/03_prediction/chemberta
./run_chemberta.sh base      # seyonec/ChemBERTa-zinc-base-v1   (~22.5 min to embed)
./run_chemberta.sh mtr       # DeepChem/ChemBERTa-77M-MTR       (~6 min to embed)
```

Each run writes the predictions (`*-Top10*.csv`), per-query metrics (`*_per.csv`) and
summary means (`Validation_*_ave.csv`) into `results_base/` or `results_mtr/`. The reported
values are in the dissertation.

> **Cache trap.** `fusion_MolTarPred.py` reuses `artifacts/cb_sim.npy` and
> `artifacts/tan_sim.npy` **without validating them**. They are `[n_query x n_db]` matrices
> keyed to the database they were built from, so after any database rebuild they must be
> **deleted**, not merely regenerated — otherwise fusion silently reads targets from
> misaligned rows.

---

## Web application

Serves MolTarPred v2.2 over the strict database (899,623 molecules, 6,758 targets,
1,349,633 ligand–target interactions). Measured: 3.1 s index load, 0.78 GB resident,
~190 ms per prediction.

```bash
# backend
cd webapp/backend
python build_index.py            # one-off, ~106 s: precompute the fingerprint index
uvicorn main:app --port 8000

# frontend (separate terminal)
cd webapp/frontend
npm install && npm run dev       # http://localhost:5173
```

**Endpoints** — `GET /api/health`, `POST /api/predict` `{smiles, top_k}`,
`GET /api/molecule` (RDKit SVG). Input is capped at 600 SMILES characters and
`top_k` ≤ 50.

The frontend returns the parsed query structure, a sortable and paginated table of
predicted targets with reliability scores and ChEMBL links, the supporting neighbours for
any selected target, and the ten nearest neighbours as structures — so a user can inspect
the evidence behind each prediction. See `webapp/README.md` for deployment notes.

---

## What is not in this repository

Roughly 14 GB of intermediate data is excluded and regenerable from the scripts above:
the ChEMBL 36 dump, the filtered and grouped databases (100 MB – 1.8 GB each), the
transformer embeddings (0.3 – 2.6 GB each), the cached similarity matrices, and the web
app's fingerprint index.

Superseded experiments are also excluded — the loose-database ChemBERTa runs, the
corrected-tokenizer sibling pipeline, and an earlier standalone embedding library. Only
the code that produced the reported results is here.

---

## Notes on the data

Three defects in the reference data pipeline were found and corrected before the reported
results were generated. They are worth knowing about if you compare numbers against
earlier MolTarPred publications:

1. **`CHEMBL612545` ("Unchecked")** is ChEMBL's placeholder for assays whose target was
   never curated — no protein components, ~2.3 M unrelated activities. It was being treated
   as a protein target, and was ground truth for 50 of the 100 benchmark drugs. Now excluded
   via `target_type != 'UNCHECKED'`.
2. **Deduplication was keyed on `compound_records.compound_key`** — the paper-local label an
   author gives a compound ("1", "5a"), shared by 47,877 distinct molecules. This collapsed
   unrelated compounds and destroyed ~497,000 real molecule–target pairs (37% of the
   interaction data). Now keyed on `molecule_chembl_id`.
3. **The SQL had no `ORDER BY`**, so deduplication kept an arbitrary row and consecutive runs
   disagreed on ~23,000 pairs. Now ordered by `activities.activity_id`; runs are
   byte-identical.

All three affected both database variants identically, and correcting them preserved every
method ranking.

---

## References

He T., Caba K., Ballester P. J. *A precise comparison of molecular target prediction
methods.* **Digital Discovery**, 2025, **4**, 2548–2558.

Peón A., Naulaerts S., Ballester P. J. *Predicting the reliability of drug–target
interaction predictions with maximum coverage of target space.* **Scientific Reports**,
2017, **7**, 3820. — the source of the Morgan + Tanimoto + *k* = 10 configuration and the
reliability score.

Peón A., Li H., Ghislat G. *et al.* *MolTarPred: A web tool for comprehensive target
prediction with reliability estimation.* **Chemical Biology & Drug Design**, 2019, **94**,
1390–1401.
