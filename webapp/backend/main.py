# main.py — FastAPI backend for MolTarPred v2.
# Run:  uvicorn main:app --reload --port 8000

from __future__ import annotations

from fastapi import FastAPI, HTTPException, Query, Response
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

import core

app = FastAPI(title="MolTarPred v2 API", version="2.0")

# CORS: allow the React dev server and same-origin production use.
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_methods=["GET", "POST"],
    allow_headers=["*"],
)

MAX_SMILES_LEN = 600  # guard against abusive input on a public endpoint


class PredictRequest(BaseModel):
    smiles: str = Field(..., min_length=1, max_length=MAX_SMILES_LEN)
    top_k: int = Field(10, ge=1, le=50)


@app.get("/api/health")
def health():
    ix = core.get_index()
    return {
        "status": "ok",
        "n_molecules": ix.n_molecules,
        "n_targets": ix.n_targets,
        "n_interactions": ix.n_interactions,
    }


@app.post("/api/predict")
def predict(req: PredictRequest):
    result = core.predict(req.smiles.strip(), top_k=req.top_k)
    if result is None:
        raise HTTPException(status_code=422, detail="Invalid SMILES string.")
    return {"query": {"smiles": req.smiles.strip(), "top_k": req.top_k}, **result}


@app.get("/api/molecule")
def molecule(
    smiles: str = Query(..., max_length=MAX_SMILES_LEN),
    w: int = Query(320, ge=80, le=900),
    h: int = Query(240, ge=80, le=900),
    legend: str = Query(""),
):
    """Render a molecule to a clean SVG (scalable, theme-friendly)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise HTTPException(status_code=422, detail="Invalid SMILES string.")
    Chem.rdDepictor.SetPreferCoordGen(True)
    Chem.rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 2
    opts.padding = 0.12
    opts.clearBackground = False           # transparent -> blends with any card colour
    opts.addStereoAnnotation = True
    if legend:
        opts.legendFontSize = 16
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol, legend=legend)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return Response(content=svg, media_type="image/svg+xml")
