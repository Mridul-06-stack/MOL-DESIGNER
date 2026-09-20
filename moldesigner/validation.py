"""Validation module – redock report, enrichment."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .scoring import Docker


@dataclass
class RedockReport:
    """Result of re-docking a molecule to verify score stability."""
    
    smiles: str
    original_scores: dict[str, float]
    redock_scores: dict[str, float]
    rmsd: float | None
    consistent: bool


def redock(
    docker: Docker,
    top_smiles: list[str],
    original_dock_scores: list[dict[str, float] | None],
    consistency_threshold: float = 2.0
) -> list[RedockReport]:
    """
    Re-dock top molecules and verify scores are consistent with originals.
    RMSD cannot be easily calculated with standard vina python API without saving PDBQTs,
    so we focus on docking energy consistency.
    """
    reports = []
    
    # Run independent re-docking
    redocked_results = docker.dock_many(top_smiles)
    
    for smi, orig, redocked in zip(top_smiles, original_dock_scores, redocked_results):
        if orig is None or redocked is None:
            continue
            
        consistent = True
        for tgt in docker.targets:
            if tgt in orig and tgt in redocked:
                if abs(orig[tgt] - redocked[tgt]) > consistency_threshold:
                    consistent = False
            else:
                consistent = False
                
        reports.append(
            RedockReport(
                smiles=smi,
                original_scores=orig,
                redock_scores=redocked,
                rmsd=None,  # RMSD calculation skipped for now
                consistent=consistent
            )
        )
        
    return reports


def enrichment_summary(reports: list[RedockReport]) -> dict[str, Any]:
    """Compute summary statistics for a batch of redock reports."""
    total = len(reports)
    if total == 0:
        return {"total": 0}
        
    consistent_count = sum(1 for r in reports if r.consistent)
    
    total_dev = 0.0
    measurements = 0
    for r in reports:
        for tgt, orig_energy in r.original_scores.items():
            if tgt in r.redock_scores:
                redock_energy = r.redock_scores[tgt]
                total_dev += abs(orig_energy - redock_energy)
                measurements += 1
                
    mean_dev = total_dev / measurements if measurements > 0 else 0.0
    
    return {
        "total": total,
        "consistent": consistent_count,
        "consistency_rate": consistent_count / total,
        "mean_deviation": mean_dev
    }
