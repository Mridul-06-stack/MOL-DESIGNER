import sys
import os
import json
import asyncio
from typing import Any
from fastapi import FastAPI
from fastapi.responses import StreamingResponse
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from rdkit import Chem

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from moldesigner.scoring import Scorer, Weights
from moldesigner.generator import evolve_islands, evolve
from moldesigner.docking import load_docker

app = FastAPI(title="MolDesigner API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class EvolveRequest(BaseModel):
    islands: int = 3
    pop_size: int = 30
    generations: int = 15
    seed: int = 42
    seed_smiles: str = "C=Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"
    targets: list[str] = ["WT", "L858R"]
    dock_weight: float = 0.5

def _json_default(obj: Any) -> Any:
    import numpy as np
    if isinstance(obj, np.floating):
        return float(obj)
    raise TypeError

@app.post("/api/evolve")
async def run_evolution(req: EvolveRequest):
    docker = load_docker({
        "backend": "rdkit_score",
        "targets": req.targets
    })
    
    scorer = Scorer(
        weights=Weights(
            dock=req.dock_weight, 
            qed=0.2, 
            sa=0.3, 
            sim=0.0
        ), 
        docker=docker
    )
    # We use a synchronous generator underlying this, so for FastAPI to stream smoothly,
    # we yield events, but wrap them in SSE format.
    def event_stream():
        if req.islands > 1:
            pop_per = max(5, req.pop_size // req.islands)
            gen = evolve_islands(
                scorer=scorer,
                seeds=[req.seed_smiles],
                islands=req.islands,
                pop_size_per_island=pop_per,
                generations=req.generations,
                seed=req.seed
            )
        else:
            gen = evolve(
                scorer=scorer,
                seeds=[req.seed_smiles],
                pop_size=req.pop_size,
                generations=req.generations,
                seed=req.seed
            )
            
        for event in gen:
            # Format as SSE
            data = json.dumps(event, default=_json_default)
            yield f"data: {data}\n\n"
            
    return StreamingResponse(event_stream(), media_type="text/event-stream")
@app.get("/api/molblock")
async def get_molblock(smiles: str):
    from rdkit.Chem import AllChem
    from fastapi.responses import PlainTextResponse
    
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return PlainTextResponse("", status_code=400)
        
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    
    block = Chem.MolToMolBlock(mol)
    return PlainTextResponse(block)

