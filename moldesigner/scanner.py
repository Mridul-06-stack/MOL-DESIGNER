"""
Escape Scanner module for adversarial red-teaming.

Takes a successful molecule and predicts a protein mutation (variant)
that would drastically reduce its binding affinity, simulating continuous
tumor evolution.
"""

from typing import Protocol, runtime_checkable
from base64 import b64encode


@runtime_checkable
class EscapeScanner(Protocol):
    """Protocol for scanning a molecule and predicting a breaking mutation."""
    
    def scan(self, smiles: str, baseline_energy: float) -> str | None:
        """
        Analyzes the SMILES and returns a new variant ID string that breaks it,
        or None if no escape mutation is predicted.
        """
        ...


class MockScanner:
    """
    Simulates an escape scanner.
    Always generates a highly specific resistance mutation against the exact
    scaffold provided, encoded into the variant name itself so the Mock Docker
    can interpret it and penalize it.
    """
    
    def scan(self, smiles: str, baseline_energy: float) -> str | None:
        """
        Generates a synthetic variant name string identifying this scaffold.
        """
        if baseline_energy > -5.0:
            # If it's already a terrible molecule, no need to mutate to escape it
            return None
            
        # Create a tiny hash representing the core scaffold/features
        # We just base64 the smiles (or first half of it) to establish a unique ID.
        # RDKitScoreDocker will look for "RESIST_SIM_" strings and use Tanimoto 
        # similarity to penalize anything resembling this candidate.
        b64 = b64encode(smiles.encode("utf-8")).decode("utf-8")
        escaped_id = b64[:8].replace("=", "").replace("+", "").replace("/", "")
        
        # We append the actual smiles behind a delimiter so the docker can parse it
        return f"RESIST_SIM_{escaped_id}_{smiles}"
