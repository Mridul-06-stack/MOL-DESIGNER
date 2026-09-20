#!/usr/bin/env python3
"""
EGFR Panel CLI for MolDesigner.

Runs a multi-variant evolutionary optimization against the EGFR clinical panel
(WT, L858R, T790M_C797S) using the EGFRMockDocker.

Demonstrates how worst-case scoring forces the evolutionary algorithm to 
seek structural classes that evade resistance (like T790M/C797S).

Usage:
    python cli/run_egfr_panel.py --generations 15 --pop-size 40
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from moldesigner.scoring import Scorer, Weights
from moldesigner.generator import evolve, evolve_islands
from moldesigner.docking import load_docker


def main():
    parser = argparse.ArgumentParser(
        description="MolDesigner EGFR Panel Run",
    )
    parser.add_argument("--pop-size", type=int, default=40, help="Population size")
    parser.add_argument("--generations", type=int, default=15, help="Number of generations")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--islands", type=int, default=1, help="Number of islands")
    args = parser.parse_args()

    # Erlotinib (a known 1st-gen EGFR inhibitor, binds WT/L858R but has resistance to T790M)
    # We will use this as our seed
    erlotinib_smiles = "C=Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"

    # Heavily weight docking so resistance forces evolution away from base scaffold
    weights = Weights(dock=0.6, qed=0.2, sa=0.2, sim=0.0)

    # Load the specialized EGFR mock docker
    docker = load_docker({"backend": "egfr_mock"})

    scorer = Scorer(
        weights=weights,
        docker=docker,
    )

    print("=== MolDesigner EGFR Resistance Panel ===")
    print(f"Targets: {docker.targets}")
    print(f"Weights: {weights.as_dict()}")
    print("Optimization focuses on the *worst-case* binding across the panel.")
    print("Seed molecule: Erlotinib (fails on T790M_C797S)")
    print(f"Islands: {args.islands}")
    print("-" * 50)

    if args.islands > 1:
        pop_per_island = max(10, args.pop_size // args.islands)
        evolution_gen = evolve_islands(
            scorer=scorer,
            seeds=[erlotinib_smiles],
            islands=args.islands,
            pop_size_per_island=pop_per_island,
            generations=args.generations,
            seed=args.seed,
        )
    else:
        evolution_gen = evolve(
            scorer=scorer,
            seeds=[erlotinib_smiles],
            pop_size=args.pop_size,
            generations=args.generations,
            seed=args.seed,
        )

    for event in evolution_gen:
        gen = event["generation"]
        best = event["best"]
        mean = event["mean_fitness"]
        new_mols = event.get("new_molecules", 0)

        if best:
            dock_worst = best["raw"]["dock_worst"]
            dock_scores = best["raw"]["dock"]
            
            # Format scores directly
            wt = dock_scores.get("WT", 0.0)
            l858r = dock_scores.get("L858R", 0.0)
            t790m = dock_scores.get("T790M_C797S", 0.0)
            
            print(f"Gen {gen:2d} | Fit: {best['fitness']:.3f} | Worst: {dock_worst:.2f}")
            print(f"       | [WT: {wt:5.2f}] [L858R: {l858r:5.2f}] [T790M: {t790m:5.2f}]")
            print(f"       | {best['smiles'][:60]}")
            print("-" * 50)


if __name__ == "__main__":
    main()
