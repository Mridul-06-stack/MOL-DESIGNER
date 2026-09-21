"""
Escape Scanner module for adversarial red-teaming.

Takes a successful molecule and predicts a protein mutation (variant)
that would drastically reduce its binding affinity, simulating continuous
tumor evolution.
"""

from dataclasses import dataclass
from typing import Protocol, runtime_checkable
from base64 import b64encode
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors


@runtime_checkable
class EscapeScanner(Protocol):
    """Protocol for scanning a molecule and predicting a breaking mutation."""
    
    def scan(self, smiles: str, baseline_energy: float) -> str | None:
        """
        Analyzes the SMILES and returns a new variant ID string that breaks it,
        or None if no escape mutation is predicted.
        """
        ...


@dataclass
class ClinicalMutation:
    id: str
    name: str
    exon: str
    mechanism: str
    clinical_context: str
    adaptation_guidance: str


CLINICAL_EGFR_MUTATIONS: dict[str, ClinicalMutation] = {
    "T790M": ClinicalMutation(
        id="T790M",
        name="Gatekeeper T790M",
        exon="Exon 20",
        mechanism="Steric Clash (+54 Å³ Met Sidechain)",
        clinical_context="Primary acquired resistance to 1st/2nd Gen TKIs (Erlotinib/Gefitinib). Found in ~50% of resistant lung adenocarcinomas.",
        adaptation_guidance="Evolve a compact core (MW < 420) or flexible linkers (sulfonamide, morpholine) that avoid the Met790 gatekeeper bump."
    ),
    "C797S": ClinicalMutation(
        id="C797S",
        name="Covalent Loss C797S",
        exon="Exon 20",
        mechanism="Abolishes Nucleophilic Thiol (Loss of Covalent Anchor)",
        clinical_context="Primary acquired resistance to 3rd Gen irreversible TKIs (Osimertinib). Cys->Ser substitution prevents covalent bond formation.",
        adaptation_guidance="Transition from covalent reliance to high-affinity reversible allosteric binding or multi-point hydrogen bond networks."
    ),
    "L718Q": ClinicalMutation(
        id="L718Q",
        name="P-Loop Hinge L718Q",
        exon="Exon 18",
        mechanism="Hydrophobic to Polar Shift (Electrostatic Repulsion)",
        clinical_context="Causes clinical cross-resistance to Osimertinib and Lazertinib in patients without T790M.",
        adaptation_guidance="Reduce hydrophobic grease (LogP < 3.2) and introduce hydrogen bond acceptors to complement the polar glutamine sidechain."
    ),
    "G724S": ClinicalMutation(
        id="G724S",
        name="ATP-Pocket G724S",
        exon="Exon 19",
        mechanism="Kinase C-Helix Structural Distortion",
        clinical_context="Conferring resistance to Osimertinib in Exon 19 deletion lung cancers.",
        adaptation_guidance="Evolve conformationally adaptive scaffolds capable of binding the strained ATP-cleft."
    ),
}


class ClinicalEGFRScanner:
    """
    Oncology-grounded escape scanner for human EGFR kinase.
    Evaluates candidate molecular features and introduces realistic clinical
    resistance mutations sequentially or based on specific chemical vulnerabilities.
    """

    def __init__(self, energy_threshold: float = -5.0) -> None:
        self.energy_threshold = energy_threshold

    def scan(self, smiles: str, baseline_energy: float, active_targets: list[str] | None = None) -> str | None:
        """
        Analyzes the lead SMILES and returns the next clinical EGFR mutation ID,
        or None if the molecule is not potent enough to trigger tumor mutation.
        """
        if baseline_energy > self.energy_threshold:
            return None

        active = set(active_targets or [])
        mol = Chem.MolFromSmiles(smiles)

        # 1. Gatekeeper T790M: primary evolutionary response to initial potent leads
        if "T790M" not in active:
            return "T790M"

        # 2. C797S: tumor response when molecules feature covalent warheads or overcome T790M
        if "C797S" not in active:
            if mol is not None:
                # Check for covalent warheads: acrylamide C=C-C(=O)N, vinyl sulfone C=C-S(=O)(=O)
                has_acrylamide = mol.HasSubstructMatch(Chem.MolFromSmarts("C=CC(=O)N"))
                has_vinyl_sulfone = mol.HasSubstructMatch(Chem.MolFromSmarts("C=CS(=O)(=O)"))
                if has_acrylamide or has_vinyl_sulfone or "T790M" in active:
                    return "C797S"
            else:
                return "C797S"

        # 3. L718Q: secondary/4th-gen resistance targeting lipophilic leads
        if "L718Q" not in active:
            return "L718Q"

        # 4. G724S: tertiary structural distortion
        if "G724S" not in active:
            return "G724S"

        return None

    def get_mutation_info(self, mutation_id: str) -> ClinicalMutation | None:
        """Retrieve rich oncological metadata for an EGFR mutation."""
        return CLINICAL_EGFR_MUTATIONS.get(mutation_id)


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
            return None
            
        b64 = b64encode(smiles.encode("utf-8")).decode("utf-8")
        escaped_id = b64[:8].replace("=", "").replace("+", "").replace("/", "")
        return f"RESIST_SIM_{escaped_id}_{smiles}"
