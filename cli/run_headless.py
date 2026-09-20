#!/usr/bin/env python3
"""
Headless CLI for MolDesigner evolution runs.

Runs the genetic algorithm without docking (Phase 1) or with docking
when configured. Outputs JSONL run logs for replay.

Usage:
    python cli/run_headless.py --seeds "CCO" "c1ccccc1" --generations 12 --seed 42
"""

from __future__ import annotations

import argparse
import json
import os
import sys

# Ensure the parent directory is on the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from moldesigner.scoring import Scorer, Weights
from moldesigner.generator import evolve, evolve_islands
from moldesigner.docking import load_docker
from moldesigner.validation import redock, enrichment_summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description="MolDesigner headless evolution run",
        epilog="Docking scores are proxies – not experimentally validated.",
    )
    parser.add_argument(
        "--seeds", nargs="+", required=True,
        help="Seed SMILES for the initial population",
    )
    parser.add_argument("--pop-size", type=int, default=60, help="Population size")
    parser.add_argument("--generations", type=int, default=30, help="Number of generations")
    parser.add_argument("--seed", type=int, default=0, help="Random seed for reproducibility")
    parser.add_argument("--islands", type=int, default=1, help="Number of island sub-populations")
    parser.add_argument(
        "--crossover-frac", type=float, default=0.5,
        help="Fraction of offspring from crossover vs mutation",
    )
    parser.add_argument(
        "--diversity-max-sim", type=float, default=0.75,
        help="Max Tanimoto similarity for diversity filter",
    )
    parser.add_argument(
        "--weights", type=str, default=None,
        help='JSON string of weights, e.g. \'{"qed":0.3,"sa":0.2}\'',
    )
    parser.add_argument(
        "--references", nargs="*", default=None,
        help="Reference SMILES for similarity scoring",
    )
    
    # Docking arguments
    parser.add_argument(
        "--dock-backend", type=str, choices=["vina", "rdkit_score"], default=None,
        help="Backend to use for docking (None disables docking)"
    )
    parser.add_argument(
        "--dock-weight", type=float, default=0.0,
        help="Weight for the docking score component"
    )
    parser.add_argument(
        "--dock-targets", nargs="+", default=["WT", "L858R"],
        help="Target variants to dock against (default: WT L858R)"
    )
    parser.add_argument(
        "--dock-config", type=str, default=None,
        help="JSON string of additional docker config (e.g., receptor_paths for vina)"
    )

    parser.add_argument(
        "--output", type=str, default=None,
        help="Path to output JSONL file (default: runs/<seed>.jsonl)",
    )

    args = parser.parse_args()

    # Build weights
    weights = Weights()
    if args.weights:
        w_dict = json.loads(args.weights)
        for k, v in w_dict.items():
            if hasattr(weights, k):
                setattr(weights, k, float(v))
                
    if args.dock_backend is not None:
        weights.dock = args.dock_weight

    # Output path
    output_path = args.output
    if output_path is None:
        os.makedirs("runs", exist_ok=True)
        output_path = os.path.join("runs", f"run_seed{args.seed}.jsonl")

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else ".", exist_ok=True)
    
    # Configure docker
    docker = None
    if args.dock_backend is not None:
        config = {"backend": args.dock_backend, "targets": args.dock_targets}
        if args.dock_config:
            config.update(json.loads(args.dock_config))
        docker = load_docker(config)

    # Build scorer
    scorer = Scorer(
        weights=weights,
        docker=docker,
        references=args.references,
    )

    print(f"MolDesigner Headless Run")
    print(f"  Seeds: {args.seeds}")
    
    if args.islands > 1:
        print(f"  Islands: {args.islands} (Pop size per island: {args.pop_size // args.islands})")
    else:
        print(f"  Pop size: {args.pop_size}")
        
    print(f"  Generations: {args.generations}")
    print(f"  Seed: {args.seed}")
    print(f"  Weights: {weights.as_dict()}")
    if docker:
        print(f"  Docker: {args.dock_backend} (Targets: {docker.targets})")
    print(f"  Output: {output_path}")
    print(f"  Disclaimer: docking scores are proxies – not experimentally validated.")
    print()

    last_event = None
    with open(output_path, "w") as f:
        if args.islands > 1:
            pop_per_island = max(10, args.pop_size // args.islands)
            evolution_gen = evolve_islands(
                scorer=scorer,
                seeds=args.seeds,
                islands=args.islands,
                pop_size_per_island=pop_per_island,
                generations=args.generations,
                crossover_frac=args.crossover_frac,
                diversity_max_sim=args.diversity_max_sim,
                seed=args.seed,
            )
        else:
            evolution_gen = evolve(
                scorer=scorer,
                seeds=args.seeds,
                pop_size=args.pop_size,
                generations=args.generations,
                crossover_frac=args.crossover_frac,
                diversity_max_sim=args.diversity_max_sim,
                seed=args.seed,
            )
            
        for event in evolution_gen:
            # Write event
            f.write(json.dumps(event, default=_json_default) + "\n")
            f.flush()

            gen = event["generation"]
            best = event["best"]
            mean = event["mean_fitness"]
            new = event.get("new_molecules", 0)
            
            last_event = event

            if best:
                print(
                    f"Gen {gen:3d} | best={best['fitness']:.4f} | "
                    f"mean={mean:.4f} | new={new} | "
                    f"{best['smiles'][:60]}"
                )
            else:
                print(f"Gen {gen:3d} | empty population")

    print(f"\nRun complete. Log saved to {output_path}")
    print(f"Cache: {scorer.cache_size} molecules scored")
    
    # Validation step
    if docker and last_event and last_event.get("top"):
        print("\n=== Validation Phase ===")
        print("Running independent re-docking of top 5 molecules...")
        top_candidates = last_event["top"][:5]
        top_smiles = [c["smiles"] for c in top_candidates]
        original_scores = [c["raw"].get("dock", {}) for c in top_candidates]
        
        reports = redock(docker, top_smiles, original_scores)
        summary = enrichment_summary(reports)
        
        print(f"Redocked {summary['total']} molecules.")
        print(f"Consistent: {summary['consistent']} ({summary['consistency_rate']*100:.1f}%)")
        print(f"Mean deviation: {summary['mean_deviation']:.2f}")


def _json_default(obj):
    """JSON serialiser fallback for numpy types."""
    import numpy as np
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.bool_):
        return bool(obj)
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


if __name__ == "__main__":
    main()
