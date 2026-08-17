# MolTarPred v2 — Web Application

Predict protein targets for a small molecule by Morgan/Tanimoto similarity search
over curated **ChEMBL 36** bioactivity data (the **MolTarPred v2** database:
binding assays, IC50/Ki/EC50 with `relation = '='` and value ≤ 10,000 nM,
excluding ChEMBL's `UNCHECKED` placeholder target).

> The v2 database was rebuilt on 2026-07-16 (Unchecked pseudo-target removed,
> deduplication corrected, SQL made deterministic) — see
> `Code/DATA_PIPELINE_FIXES.md`. It grew from 680k to ~900k molecules and from
> 964k to ~1.35M ligand-target interactions. **Re-run `build_index.py` after any
> such rebuild** or the app silently serves the old data.

The app is split into a **FastAPI backend** (prediction API) and a **React SPA**
frontend. The original Streamlit app (`app.py`, `moltarpred.py`) is retained in
this folder for reference.

```
webapp/fastapi_react/
├── backend/            FastAPI service
│   ├── main.py         API routes (predict, molecule SVG, health)
│   ├── core.py         framework-agnostic prediction logic (v2)
│   ├── build_index.py  one-off: precompute fingerprints + target map
│   ├── requirements.txt
│   └── data/           generated artifacts (git-ignored, ~150 MB on disk)
└── frontend/           Vite + React SPA
    ├── src/            App.jsx, api.js, components/, styles.css
    └── package.json
```

> The legacy Streamlit app lives alongside this one in `webapp/streamlit/`.

## Architecture

- **Backend** loads a precomputed fingerprint index once at startup (`core.Index`,
  memoised) so there is no ~105 s fingerprinting cost per boot. Measured: 3.1 s to
  load, 0.78 GB resident, ~190 ms per prediction. Endpoints:
  - `GET  /api/health` → molecule / target / interaction counts
  - `POST /api/predict` `{smiles, top_k}` → ranked targets (+ supporting ligands) & top hits
  - `GET  /api/molecule?smiles=&w=&h=&legend=` → clean transparent-background SVG
- **Frontend** calls same-origin `/api/*`; in dev, Vite proxies to the backend
  (`vite.config.js`).

## Prerequisites

- Python env with `rdkit`, `pandas` (the project's `chem-gnn-gen` conda env) plus
  `fastapi` + `uvicorn` (`pip install -r backend/requirements.txt`).
- Node 18+ / npm.

## Run (development)

**1. Build the search index once** (reads the v2 CSVs, writes `backend/data/`):

```bash
cd backend
python build_index.py          # ~105 s, fingerprints ~900k molecules
```

**2. Start the backend:**

```bash
cd backend
uvicorn main:app --port 8000   # add --reload while developing
```

**3. Start the frontend:**

```bash
cd frontend
npm install
npm run dev                     # http://localhost:5173
```

Open http://localhost:5173.

## Production build

```bash
cd frontend && npm run build    # emits static assets to frontend/dist/
```

Serve `frontend/dist/` from any static host (or mount it in FastAPI with
`StaticFiles`) and run the backend behind a process manager. For public traffic:

- Run **N uvicorn workers / replicas** behind a reverse proxy (nginx / Caddy);
  each replica holds the index in RAM — **~0.78 GB resident** (measured), not the
  ~150 MB on-disk pickle size. Size accordingly.
- The `/api/predict` endpoint already caps SMILES length and `top_k`
  (`main.py`); add rate limiting at the proxy.
- Regenerate `backend/data/` by re-running `build_index.py` whenever the v2
  database changes.

## Updating the database

`build_index.py` points at:
- `…/Code/v2/v2_chembl_36_database.csv` (molecules) — the **unmasked** database
- `…/Code/v2/v2_chembl_36_filter.csv` (target names / organisms)

Edit those paths to swap in a different MolTarPred database, then re-run it.

### Masked vs unmasked — do not "fix" this

The served index is built from the **unmasked** database and must stay that way.
`v2_chembl_36_database_masked.csv` has the 100 benchmark query ligands removed so
that a query cannot retrieve itself during validation. That is required for the
reported metrics and wrong for the deployed tool: those 100 are approved drugs,
and withholding them silently impoverishes the neighbourhood returned to any user
searching near that chemistry — including hiding 5 targets that no other molecule
in the database carries.

Masking belongs to the evaluation pipeline only. The consequence of serving the
unmasked database is that querying a drug already in ChEMBL returns that drug as
its own top hit at similarity 1.0, which is correct behaviour for a lookup tool.
