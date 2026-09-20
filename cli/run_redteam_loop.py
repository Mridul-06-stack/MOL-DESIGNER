#!/usr/bin/env python3
"""
Phase 5 CLI — Red-Team Escape Loop

Executes a cyclical pipeline:
1. Evolve molecules against current targets (starting with WT).
2. Take the best candidate and ask the EscapeScanner for a resistance mutation.
3. Add the mutation to the target list.
4. Repeat. Watch the generator get pushed to explore new chemistry!
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from moldesigner.scoring import Scorer, Weights
from moldesigner.generator import evolve_islands
from moldesigner.docking import RDKitScoreDocker
from moldesigner.scanner import MockScanner


def main():
    print("=== MolDesigner Red-Team Loop ===")
    
    seeds = ["C=Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"] # Erlotinib
    
    # Initialize components
    scanner = MockScanner()
    targets = ["WT", "L858R"]
    
    cycles = 3
    generations_per_cycle = 5
    pop_size = 30
    
    for cycle in range(1, cycles + 1):
        print(f"\n--- CYCLE {cycle} ---")
        print(f"Current Targets: {targets}")
        
        # We supply the *updated* target list to the Docker
        docker = RDKitScoreDocker(targets=targets)
        scorer = Scorer(weights=Weights(dock=0.5, qed=0.2, sa=0.3, sim=0.0), docker=docker)
        
        # Run island evolution
        print(f"Evolving for {generations_per_cycle} generations...")
        evolution_gen = evolve_islands(
            scorer=scorer,
            seeds=seeds,
            islands=3,
            pop_size_per_island=pop_size // 3,
            generations=generations_per_cycle,
            seed=42 + cycle
        )
        
        final_gen = None
        for step in evolution_gen:
            final_gen = step
            
        best = final_gen["best"]
        print(colored_best(best, targets, cycle))
        
        # Run Red-Team Scanner Phase
        print("\n[Scanner Phase]")
        # We pass the best molecule's SMILES and worst-case energy to the scanner.
        worst_energy = max(best["parts"].get(f"raw_{t}", 0.0) for t in targets)
        
        new_variant = scanner.scan(best["smiles"], worst_energy)
        if new_variant:
            print(f"🚨 SCANNER ALERT! Mutational vulnerability discovered.")
            print(f"  New Target Injected: {new_variant[:30]}...")
            
            # The new variant gets permanently appended to the board
            targets.append(new_variant)
            
            # We seed the next cycle with the best molecules we had 
            # to watch them instantly fail and be forced to evolve
            seeds = [s["smiles"] for s in final_gen["top"][:3]]
        else:
            print(f"✅ Scanner could not find a plausible escape mutation. Scaffold robust.")
            break


def colored_best(best, targets, cycle):
    smiles = best['smiles']
    fit = best['fitness']
    
    stats = []
    for tgt in targets:
        # raw value is exposed dynamically as raw_WT, etc.
        val = best["parts"].get(f"raw_{tgt}", 0.0)
        stats.append(f"[{tgt[:15]}: {val:.2f}]")
        
    worst = max(best["parts"].get(f"raw_{tgt}", 0.0) for t in targets)
    return (
        f"Result Gen {cycle * 5} | Fit: {fit:.3f} | Worst: {worst:.2f}\n"
        f"       | {' '.join(stats)}\n"
        f"       | {smiles}"
    )

if __name__ == "__main__":
    main()
