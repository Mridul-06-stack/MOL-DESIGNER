"""Tests for moldesigner.generator evolve_islands."""

import pytest
from rdkit import Chem

from moldesigner.scoring import Scorer, Weights
from moldesigner.generator import evolve_islands

ERLOTINIB = "C=Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"

def test_evolve_islands_runs():
    """Verify that the island model successfully iterates through generations."""
    scorer = Scorer(weights=Weights())
    
    gen = evolve_islands(
        scorer=scorer,
        seeds=[ERLOTINIB],
        islands=2,
        pop_size_per_island=5,
        generations=2,
        migration_interval=1,
        seed=42
    )
    
    events = list(gen)
    assert len(events) == 3  # Gen 0, 1, 2
    
    for event in events:
        assert "generation" in event
        assert "best" in event
        assert "top" in event
        assert len(event["top"]) > 0


def test_islands_merge():
    """Verify that multiple islands merge their results correctly at the end."""
    scorer = Scorer(weights=Weights(qed=1.0))
    gen = evolve_islands(
        scorer=scorer,
        seeds=[ERLOTINIB, ASPIRIN],
        islands=3,
        pop_size_per_island=10,
        generations=2,
        migration_interval=2, 
        seed=123
    )

    
    last_event = list(gen)[-1]
    assert last_event["best"] is not None
    assert last_event["mean_fitness"] > 0
