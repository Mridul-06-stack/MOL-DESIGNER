"""
LLM Agent Interface for MolDesigner.

This module exposes functional paradigms designed to be ingested by
tool-calling LLMs (OpenAI, Claude, LangChain, etc.). 
Instead of operating the CLI, an AI can autonomously trigger 
evolution targets and evaluate the red-teaming outputs.
"""

from typing import Dict, Any, List
import json
import asyncio

from moldesigner.scoring import Scorer, Weights
from moldesigner.generator import evolve_islands
from moldesigner.scanner import MockScanner
from moldesigner.docking import load_docker

async def run_evolution_task(
    target_names: List[str], 
    aggressiveness: float = 0.5, 
    seed_molecule: str = "C=Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"
) -> Dict[str, Any]:
    """
    Agent tool designed for an LLM to call when a user says something like:
    "Evolve a drug for EGFR targeting WT and L858R. Make the search broad."
    
    Args:
        target_names: The protein variants (e.g. ["WT", "L858R", "T790M_C797S"])
        aggressiveness: A 0.0 - 1.0 float abstracting search diversity. 
                        Higher means more islands, larger populations, longer search.
        seed_molecule: Initial molecular scaffold SMILES.
    
    Returns:
        JSON-serializable dict of the final best candidate and the red-team analysis.
    """
    # LLM param abstraction mappings
    islands = max(1, int(aggressiveness * 5))
    generations = max(5, int(aggressiveness * 25))
    pop_size = max(10, int(aggressiveness * 60))
    
    docker = load_docker({"backend": "rdkit_score", "targets": target_names})
    scorer = Scorer(weights=Weights(dock=0.6, qed=0.2, sa=0.2), docker=docker)
    
    gen = evolve_islands(
        scorer, 
        [seed_molecule], 
        islands=islands, 
        pop_size_per_island=pop_size // islands, 
        generations=generations
    )
    
    final_event = None
    for event in gen:
        final_event = event
        
    best = final_event["best"]
    
    # Run the automated red-team scanner silently to return to the agent
    worst_energy = max(best["parts"].get(f"raw_{t}", 0.0) for t in target_names)
    scanner = MockScanner()
    escape_target = scanner.scan(best["smiles"], worst_energy)
    
    return {
        "status": "success",
        "best_candidate_smiles": best["smiles"],
        "simulated_affinity": worst_energy,
        "escape_mutation_vulnerability": escape_target if escape_target else "None discovered",
        "llm_advice": (
            "If an escape mutation was discovered, recursively call this tool "
            "again appending the new vulnerability string to `target_names`."
        )
    }
