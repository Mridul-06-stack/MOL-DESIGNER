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
from moldesigner.scanner import ClinicalEGFRScanner, MockScanner

app = FastAPI(title="Alchemist API")

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

_molblock_cache: dict[str, str] = {}

def _json_default(obj: Any) -> Any:
    import numpy as np
    if isinstance(obj, (np.floating, float)):
        return float(obj)
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    raise TypeError

@app.post("/api/evolve")
async def run_evolution(req: EvolveRequest):
    # Safety bounds for reliable, responsive execution on cloud tiers
    safe_islands = min(max(req.islands, 1), 3)
    safe_pop = min(max(req.pop_size, 10), 25)
    safe_gens = min(max(req.generations, 1), 30)

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

    seeds = [req.seed_smiles]
    scanner = ClinicalEGFRScanner(energy_threshold=req.redteam_threshold) if req.enable_redteam else None

    async def event_stream():
        loop = asyncio.get_running_loop()
        
        # Instantiate generator object
        gen = evolve_islands(
            scorer=scorer,
            seeds=seeds,
            islands=safe_islands,
            pop_size_per_island=safe_pop,
            generations=safe_gens,
            seed=req.seed
        )

        def get_next(iterator):
            try:
                return next(iterator), False
            except StopIteration:
                return None, True
            except Exception as e:
                return {"error": str(e)}, True

        mutations_injected = 0
        last_mutation_gen = 0

        while True:
            # Run CPU-bound generation computation in executor thread so asyncio loop handles health checks
            event, done = await loop.run_in_executor(None, get_next, gen)
            if done or event is None or "error" in event:
                break

            # Attach active targets to generation event
            event["active_targets"] = list(docker.targets)
            
            # Format as SSE
            data = json.dumps(event, default=_json_default)
            yield f"data: {data}\n\n"
            
            # Give event loop time to handle ping/health checks & network buffer flush
            await asyncio.sleep(0.03)

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
                        new_variant = scanner.scan(
                            best_mol["smiles"], 
                            dock_worst, 
                            active_targets=list(docker.targets)
                        )
                        
                        if new_variant and new_variant not in docker.targets:
                            mutations_injected += 1
                            last_mutation_gen = gen_idx
                            if hasattr(docker, "add_target"):
                                docker.add_target(new_variant)
                            else:
                                docker._targets.append(new_variant)

                            mut_info = scanner.get_mutation_info(new_variant)
                            mut_name = mut_info.name if mut_info else new_variant
                            mechanism = mut_info.mechanism if mut_info else "Biophysical steric / charge resistance"
                            context = mut_info.clinical_context if mut_info else "Clinical acquired resistance"
                            guidance = mut_info.adaptation_guidance if mut_info else "Adapt core scaffold"
                            exon = mut_info.exon if mut_info else "Exon 20"
                            
                            alert_event = {
                                "type": "redteam_alert",
                                "generation": gen_idx,
                                "mutation": new_variant,
                                "mutation_short": new_variant,
                                "mutation_name": mut_name,
                                "exon": exon,
                                "mechanism": mechanism,
                                "clinical_context": context,
                                "adaptation_guidance": guidance,
                                "trigger_smiles": best_mol["smiles"],
                                "affinity_before": dock_worst,
                                "active_targets": list(docker.targets),
                                "message": f"Clinical resistance escape: {mut_name} ({exon}). {mechanism}."
                            }
                            alert_data = json.dumps(alert_event, default=_json_default)
                            yield f"data: {alert_data}\n\n"
                            await asyncio.sleep(0.02)
            
    return StreamingResponse(event_stream(), media_type="text/event-stream")

@app.get("/api/molblock")
async def get_molblock(smiles: str):
    from rdkit.Chem import AllChem
    from fastapi.responses import PlainTextResponse
    
    if not smiles:
        return PlainTextResponse("", status_code=400)

    if smiles in _molblock_cache:
        return PlainTextResponse(_molblock_cache[smiles])

    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return PlainTextResponse("", status_code=400)
            
        mol_h = Chem.AddHs(mol)
        embed_status = AllChem.EmbedMolecule(mol_h, randomSeed=42)
        if embed_status != 0:
            embed_status = AllChem.EmbedMolecule(mol_h, useRandomCoords=True, randomSeed=42)
            
        if embed_status == 0:
            try:
                AllChem.MMFFOptimizeMolecule(mol_h, maxIters=150)
            except Exception:
                pass
            block = Chem.MolToMolBlock(mol_h)
        else:
            # Fallback to 2D
            AllChem.Compute2DCoords(mol)
            block = Chem.MolToMolBlock(mol)
            
        if len(_molblock_cache) > 300:
            _molblock_cache.clear()
        _molblock_cache[smiles] = block
        return PlainTextResponse(block)
    except Exception:
        try:
            m = Chem.MolFromSmiles(smiles)
            if m:
                AllChem.Compute2DCoords(m)
                return PlainTextResponse(Chem.MolToMolBlock(m))
        except Exception:
            pass
        return PlainTextResponse("", status_code=400)


@app.get("/api/health")
async def health_check():
    return {"status": "ok", "version": "0.1.0"}


# Serve static frontend if ui/dist exists (production/Docker unified deployment)
import os
from fastapi.staticfiles import StaticFiles

ui_dist = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "ui", "dist")
if os.path.isdir(ui_dist):
    app.mount("/", StaticFiles(directory=ui_dist, html=True), name="frontend")


