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
from moldesigner.scanner import MockScanner

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
    enable_redteam: bool = False
    redteam_threshold: float = -8.0
    redteam_interval: int = 4
    redteam_max_mutations: int = 2

def _json_default(obj: Any) -> Any:
    import numpy as np
    if isinstance(obj, np.floating):
        return float(obj)
    raise TypeError

@app.post("/api/evolve")
async def run_evolution(req: EvolveRequest):
    docker = load_docker({
        "backend": "rdkit_score",
        "targets": list(req.targets)
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

    scanner = MockScanner() if req.enable_redteam else None

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
            
        mutations_injected = 0
        last_mutation_gen = 0

        for event in gen:
            # Attach active targets to generation event
            event["active_targets"] = list(docker.targets)
            
            # Format as SSE
            data = json.dumps(event, default=_json_default)
            yield f"data: {data}\n\n"

            # Check if Red-Team should scan this generation
            if req.enable_redteam and scanner and event.get("best"):
                gen_idx = event.get("generation", 0)
                best_mol = event["best"]
                
                # Check criteria: minimum warmup, interval, cap, and good binding affinity
                if (
                    gen_idx >= 2
                    and (gen_idx - last_mutation_gen) >= req.redteam_interval
                    and mutations_injected < req.redteam_max_mutations
                ):
                    dock_worst = best_mol.get("raw", {}).get("dock_worst")
                    if dock_worst is not None and dock_worst <= req.redteam_threshold:
                        new_variant = scanner.scan(best_mol["smiles"], dock_worst)
                        
                        if new_variant and new_variant not in docker.targets:
                            mutations_injected += 1
                            last_mutation_gen = gen_idx
                            if hasattr(docker, "add_target"):
                                docker.add_target(new_variant)
                            else:
                                docker._targets.append(new_variant)

                            # Extract short clean name for UI
                            parts = new_variant.split("_", 3)
                            short_name = f"RESIST_{parts[2]}" if len(parts) >= 3 else new_variant[:15]
                            
                            alert_event = {
                                "type": "redteam_alert",
                                "generation": gen_idx,
                                "mutation": new_variant,
                                "mutation_short": short_name,
                                "trigger_smiles": best_mol["smiles"],
                                "affinity_before": dock_worst,
                                "active_targets": list(docker.targets),
                                "message": f"Resistance loophole detected! Tumor adapted with escape mutation {short_name}."
                            }
                            alert_data = json.dumps(alert_event, default=_json_default)
                            yield f"data: {alert_data}\n\n"
            
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


@app.get("/api/health")
async def health_check():
    return {"status": "ok", "version": "0.1.0"}

