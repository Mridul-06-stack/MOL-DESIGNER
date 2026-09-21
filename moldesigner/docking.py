"""Docking module – Meeko prep, Vina panel, redocking."""

from __future__ import annotations

import random
from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem

from .scoring import Docker, properties


class MeekoPrep:
    """Ligand preparation helper for AutoDock Vina via Meeko."""

    @staticmethod
    def prepare(smiles: str) -> str | None:
        """SMILES -> 3D conformer -> PDBQT string. Returns None on failure."""
        try:
            from meeko import MoleculePreparation
        except ImportError as e:
            raise ImportError(
                "Meeko is required for preparing ligands for Vina. "
                "Install it with `pip install meeko`"
            ) from e

        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
                
            # Add hydrogens and generate 3D coordinates
            mol = Chem.AddHs(mol)
            res = AllChem.EmbedMolecule(mol, randomSeed=42)
            if res != 0:
                # Fallback to random coordinates if standard embedding fails
                res = AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
                if res != 0:
                    return None
                    
            # Optimize geometry
            try:
                AllChem.MMFFOptimizeMolecule(mol)
            except Exception:
                pass  # Use unoptimized coordinates if MMFF fails
                
            preparator = MoleculePreparation(keep_nonpolar_hydrogens=False)
            preparator.prepare(mol)
            pdbqt_string = preparator.write_pdbqt_string()
            return pdbqt_string
        except Exception:
            return None


class VinaDocker:
    """Full AutoDock Vina integration conforming to the Docker protocol."""

    def __init__(
        self,
        receptor_pdbqt_paths: dict[str, str],
        center: tuple[float, float, float],
        box_size: tuple[float, float, float] = (20.0, 20.0, 20.0),
        exhaustiveness: int = 8,
        n_poses: int = 1,
    ):
        try:
            import vina
        except ImportError as e:
            raise ImportError(
                "AutoDock Vina python bindings are required. "
                "Install with `pip install vina` (may require conda-forge)"
            ) from e
            
        self.receptor_paths = receptor_pdbqt_paths
        self._targets = list(receptor_pdbqt_paths.keys())
        self.center = center
        self.box_size = box_size
        self.exhaustiveness = exhaustiveness
        self.n_poses = n_poses

    @property
    def targets(self) -> list[str]:
        return self._targets

    def dock_many(self, smiles_list: list[str]) -> list[dict[str, float] | None]:
        """Dock a batch of molecules. Returns per-mol dict[variant->energy] or None."""
        import vina
        
        results: list[dict[str, float] | None] = []
        for smi in smiles_list:
            if not smi:
                results.append(None)
                continue
                
            pdbqt = MeekoPrep.prepare(smi)
            if not pdbqt:
                results.append(None)
                continue
                
            mol_energies: dict[str, float] = {}
            success = True
            
            for variant, receptor_path in self.receptor_paths.items():
                try:
                    v = vina.Vina(sf_name='vina')
                    v.set_receptor(receptor_path)
                    v.set_ligand_from_string(pdbqt)
                    v.compute_vina_maps(center=self.center, box_size=self.box_size)
                    v.dock(exhaustiveness=self.exhaustiveness, n_poses=self.n_poses)
                    energies = v.energies(n_poses=1)
                    if energies and len(energies) > 0:
                        mol_energies[variant] = float(energies[0][0])
                    else:
                        success = False
                        break
                except Exception:
                    success = False
                    break
                    
            if success and len(mol_energies) == len(self.targets):
                results.append(mol_energies)
            else:
                results.append(None)
                
        return results


class RDKitScoreDocker:
    """
    Lightweight fallback docker using RDKit properties as proxy docking scores.
    Lets us test the pipeline end-to-end without installing Vina.
    """

    def __init__(self, targets: list[str] = None):
        self._targets = targets or ["WT", "L858R"]
        self._rng = random.Random(42)

    @property
    def targets(self) -> list[str]:
        return self._targets

    def add_target(self, target: str) -> None:
        """Dynamically register a new variant into the docking target panel."""
        if target not in self._targets:
            self._targets.append(target)

    def dock_many(self, smiles_list: list[str]) -> list[dict[str, float] | None]:
        results: list[dict[str, float] | None] = []
        for smi in smiles_list:
            try:
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    results.append(None)
                    continue
                    
                props = properties(mol)
                
                # Synthetic energy base score (more drug-like = better pseudo-energy)
                # QED is 0..1, SA is 1..10
                # Higher QED -> lower (better) energy
                # Lower SA -> lower (better) energy
                qed_term = props.get("qed", 0.5) * 5.0
                sa_term = ((10.0 - props.get("sa", 5.0)) / 9.0) * 5.0
                base_energy = -(qed_term + sa_term) # Range approx -10 to 0
                
                has_acrylamide = mol.HasSubstructMatch(Chem.MolFromSmarts("C=CC(=O)N"))
                has_vinyl_sulfone = mol.HasSubstructMatch(Chem.MolFromSmarts("C=CS(=O)(=O)"))
                is_covalent = has_acrylamide or has_vinyl_sulfone
                
                has_sulfonamide = mol.HasSubstructMatch(Chem.MolFromSmarts("S(=O)(=O)N"))
                has_morpholine = mol.HasSubstructMatch(Chem.MolFromSmarts("N1CCOCC1"))
                has_piperazine = mol.HasSubstructMatch(Chem.MolFromSmarts("N1CCNCC1"))
                has_t790m_breaker = has_sulfonamide or has_morpholine or has_piperazine or props.get("mw", 500) < 410

                mol_energies = {}
                for tgt in self.targets:
                    seed = hash(f"{smi}_{tgt}")
                    noise = random.Random(seed).uniform(-0.3, 0.3)

                    # 1. Backwards-compatible synthetic Red-Team target
                    if tgt.startswith("RESIST_SIM_"):
                        parts = tgt.split("_", 3)
                        if len(parts) >= 4:
                            target_smi = parts[3]
                            target_mol = Chem.MolFromSmiles(target_smi)
                            if target_mol is not None:
                                from rdkit.Chem import rdFingerprintGenerator
                                from rdkit import DataStructs
                                mfp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
                                fp_target = mfp_gen.GetFingerprint(target_mol)
                                fp_mol = mfp_gen.GetFingerprint(mol)
                                sim = DataStructs.TanimotoSimilarity(fp_target, fp_mol)
                                if sim > 0.6:
                                    mol_energies[tgt] = -3.0 + random.Random(seed).uniform(-0.5, 0.5)
                                    continue

                    # 2. Clinical Variant Biophysical Modeling
                    if tgt == "WT":
                        # Covalent bonus if warhead can react with Cys797
                        cov_bonus = -1.0 if is_covalent else 0.0
                        mol_energies[tgt] = base_energy + cov_bonus + noise

                    elif tgt == "L858R":
                        # Activating oncogenic driver - active site is thermodynamically primed
                        cov_bonus = -1.2 if is_covalent else 0.0
                        mol_energies[tgt] = base_energy - 0.4 + cov_bonus + noise

                    elif tgt == "T790M":
                        # Gatekeeper mutation: Met790 introduces a +54 Å³ steric clash.
                        # Compact cores or flexible hinge-binders bypass clash.
                        cov_bonus = -1.0 if is_covalent else 0.0
                        if has_t790m_breaker:
                            mol_energies[tgt] = base_energy - 0.1 + cov_bonus + noise
                        else:
                            # Bulky rigid scaffold clashes with Met790 sidechain
                            mol_energies[tgt] = base_energy + 3.2 + cov_bonus + noise

                    elif tgt == "C797S":
                        # Covalent loss: Cys797 -> Ser797 abolishes the reactive nucleophile.
                        if is_covalent:
                            mol_energies[tgt] = base_energy + 2.5 + noise
                        else:
                            mol_energies[tgt] = base_energy + 0.2 + noise

                    elif tgt == "L718Q":
                        # P-loop hydrophobic to polar shift: penalizes lipophilic grease
                        logp = props.get("logp", 3.0)
                        if logp > 3.2:
                            mol_energies[tgt] = base_energy + (logp - 3.2) * 1.5 + noise
                        else:
                            mol_energies[tgt] = base_energy - 0.2 + noise

                    elif tgt == "G724S":
                        # ATP pocket structural cleft deformation
                        rotb = props.get("rotb", 4)
                        if rotb < 3:
                            mol_energies[tgt] = base_energy + 1.8 + noise
                        else:
                            mol_energies[tgt] = base_energy + noise

                    else:
                        mol_energies[tgt] = base_energy + noise
                    
                results.append(mol_energies)
            except Exception:
                results.append(None)
                
        return results


class EGFRMockDocker:
    """
    Specialised mock docker for the EGFR panel use-case: WT, L858R, T790M_C797S.
    Simulates clinical resistance profiles to enable end-to-end multi-variant testing.
    """

    def __init__(self) -> None:
        self._targets = ["WT", "L858R", "T790M_C797S"]

    @property
    def targets(self) -> list[str]:
        return self._targets

    def dock_many(self, smiles_list: list[str]) -> list[dict[str, float] | None]:
        results: list[dict[str, float] | None] = []
        for smi in smiles_list:
            try:
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    results.append(None)
                    continue

                props = properties(mol)
                
                # Base energy driven by drug-likeness
                qed_score = props.get("qed", 0.0)
                # Cap base energy at -9.5
                base = -5.0 - (qed_score * 4.5) 
                
                # Introduce deterministic noise
                rng = random.Random(hash(smi))

                energies = {}
                
                # WT and L858R: usually sensitive to standard scaffolds
                energies["WT"] = base + rng.uniform(-0.5, 0.5)
                energies["L858R"] = base - 0.2 + rng.uniform(-0.5, 0.5)  # Slightly more sensitive

                # T790M_C797S double mutant: highly resistant
                # We simulate escaping this resistance if the molecule has an amide or specific group.
                # If it doesn't have an amide or sulfonamide, penalise heavily (steric clash / lost covalent).
                amide_smarts = Chem.MolFromSmarts("C(=O)N")
                sulfonamide = Chem.MolFromSmarts("S(=O)(=O)N")
                
                has_breaker = False
                if mol.HasSubstructMatch(amide_smarts) or mol.HasSubstructMatch(sulfonamide):
                    # Also needs to be small enough to fit the mutated pocket
                    if props.get("mw", 500) < 450:
                        has_breaker = True

                if has_breaker:
                    # Breaker scaffold: good binding
                    energies["T790M_C797S"] = base + 0.5 + rng.uniform(-0.5, 0.5)
                else:
                    # Standard scaffold: terrible binding (e.g. -4.0 kcal/mol)
                    energies["T790M_C797S"] = -4.0 + rng.uniform(-0.2, 0.2)
                    
                results.append(energies)
            except Exception:
                results.append(None)
                
        return results


def load_docker(config: dict[str, Any]) -> Docker:
    """Factory to instantiate the appropriate Docker."""
    backend = config.get("backend", "rdkit_score")
    
    if backend == "rdkit_score":
        targets = config.get("targets", ["WT", "L858R"])
        return RDKitScoreDocker(targets=targets)
        
    elif backend == "egfr_mock":
        return EGFRMockDocker()
        
    elif backend == "vina":
        receptor_paths = config.get("receptor_paths", {})
        if not receptor_paths:
            raise ValueError("VinaDocker requires receptor_paths configuration")
        center = config.get("center", (0.0, 0.0, 0.0))
        box_size = config.get("box_size", (20.0, 20.0, 20.0))
        return VinaDocker(
            receptor_pdbqt_paths=receptor_paths,
            center=center,
            box_size=box_size,
        )
        
    elif backend == "gnina":
        return GNINADocker(config.get("receptor_paths", {}))
        
    else:
        raise ValueError(f"Unknown docking backend: {backend}")

class GNINADocker:
    """
    Placeholder: Phase 7 Stretch Goal
    Translates docking requests to a local GNINA binary or external GPU server 
    to evaluate CNN-based absolute binding affinities and improved scoring.
    """
    def __init__(self, receptor_paths: dict[str, str]):
        self._targets = list(receptor_paths.keys())

    @property
    def targets(self) -> list[str]:
        return self._targets

    def dock_many(self, smiles_list: list[str]) -> list[dict[str, float] | None]:
        raise NotImplementedError("GNINA runtime integration requires Nvidia CUDA + Heavy Binary. Hook structurally prepared.")

class AiZynthFilter:
    """
    Placeholder: Phase 7 Stretch Goal
    Integrates with AiZynthFinder to parse proposed structures and evaluate
    if they are practically synthesizable in a wet lab using stock templates.
    """
    def is_synthesizable(self, smiles: str) -> bool:
        raise NotImplementedError("AiZynthFinder requires tree-search models. Hook structurally prepared.")
