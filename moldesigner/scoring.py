"""
Molecular scoring: properties, weights, fitness, and caching.

Docking scores are proxies – never label them "binding affinity".
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from typing import Any, NamedTuple, Protocol, runtime_checkable

import numpy as np
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors, QED
from rdkit.Chem import RDConfig
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

# ── SA scorer (RDKit contrib) ──────────────────────────────────────────
_sa_dir = os.path.join(RDConfig.RDContribDir, "SA_Score")
if _sa_dir not in sys.path:
    sys.path.insert(0, _sa_dir)
import sascorer  # type: ignore  # noqa: E402

# ── PAINS filter ───────────────────────────────────────────────────────
_pains_params = FilterCatalogParams()
_pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
_PAINS_CATALOG = FilterCatalog(_pains_params)


# ── Morgan fingerprint for similarity ─────────────────────────────────
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs

_mfp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def _morgan_fp(mol: Chem.Mol):
    """Morgan fingerprint (radius 2, 2048 bits)."""
    return _mfp_gen.GetFingerprint(mol)


# ── Public API ─────────────────────────────────────────────────────────

def properties(mol: Chem.Mol) -> dict[str, Any]:
    """
    Compute drug-likeness properties for a molecule.

    Returns dict with keys: qed, sa, mw, logp, hbd, hba, tpsa, rotb,
    lipinski_violations, pains.
    """
    mw = Descriptors.ExactMolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    qed_val = QED.qed(mol)
    sa_val = sascorer.calculateScore(mol)
    pains = _PAINS_CATALOG.HasMatch(mol)

    # Lipinski violations: MW>500, logP>5, HBD>5, HBA>10
    violations = sum([
        mw > 500,
        logp > 5,
        hbd > 5,
        hba > 10,
    ])

    return {
        "qed": qed_val,
        "sa": sa_val,
        "mw": mw,
        "logp": logp,
        "hbd": hbd,
        "hba": hba,
        "tpsa": tpsa,
        "rotb": rotb,
        "lipinski_violations": violations,
        "pains": pains,
    }


@dataclass
class Weights:
    """Mutable scoring weights. Each in [0, 1]."""

    dock: float = 0.0
    qed: float = 0.25
    sa: float = 0.15
    sim: float = 0.0

    def as_dict(self) -> dict[str, float]:
        return {"dock": self.dock, "qed": self.qed, "sa": self.sa, "sim": self.sim}


class Scored(NamedTuple):
    """Result of scoring a single molecule."""

    smiles: str
    fitness: float
    parts: dict[str, float]       # normalised component scores
    raw: dict[str, Any]           # raw properties + dock dict


@runtime_checkable
class Docker(Protocol):
    """Protocol for docking backends."""

    targets: list[str]

    def dock_many(self, smiles: list[str]) -> list[dict[str, float] | None]:
        """Dock a batch; returns per-mol dict[variant→energy] or None on failure."""
        ...


class Scorer:
    """
    Multi-objective scorer with aggressive caching.

    Caches *raw* components per canonical SMILES. Fitness is recomputed
    from cached components on every call so that changing weights never
    triggers re-docking.
    """

    def __init__(
        self,
        weights: Weights,
        docker: Docker | None = None,
        references: list[str] | None = None,
        dock_prefilter_max_violations: int = 1,
    ):
        self.weights = weights
        self.docker = docker
        self.dock_prefilter_max_violations = dock_prefilter_max_violations

        # Pre-compute reference fingerprints for similarity
        self._ref_fps: list = []
        if references:
            for smi in references:
                m = Chem.MolFromSmiles(smi)
                if m is not None:
                    self._ref_fps.append(_morgan_fp(m))

        # Cache: canonical SMILES → {"props": dict, "sim": float, "dock": dict|None}
        self._cache: dict[str, dict[str, Any]] = {}

    # ── Cache helpers ──────────────────────────────────────────────────

    def _get_cached(self, smi: str) -> dict[str, Any] | None:
        return self._cache.get(smi)

    def _ensure_cached(self, mol: Chem.Mol, smi: str) -> dict[str, Any]:
        """Compute and cache raw components if not already cached."""
        cached = self._cache.get(smi)
        if cached is not None:
            return cached

        props = properties(mol)

        # Similarity
        sim = 0.0
        if self._ref_fps:
            fp = _morgan_fp(mol)
            sim = max(DataStructs.TanimotoSimilarity(fp, ref) for ref in self._ref_fps)

        entry: dict[str, Any] = {"props": props, "sim": sim, "dock": None}
        self._cache[smi] = entry
        return entry

    def clear_dock_cache(self) -> None:
        """Clear cached docking results when the target panel changes."""
        for entry in self._cache.values():
            entry["dock"] = None

    # ── Batch docking ──────────────────────────────────────────────────

    def _dock_batch(self, mols_and_smiles: list[tuple[Chem.Mol, str]]) -> None:
        """Dock molecules that need docking and haven't been docked yet."""
        if self.docker is None or self.weights.dock <= 0:
            return

        to_dock: list[str] = []
        to_dock_smiles: list[str] = []
        for _mol, smi in mols_and_smiles:
            entry = self._cache.get(smi)
            if entry is None:
                continue
            # Skip if already docked against all current targets
            if entry["dock"] is not None and (
                self.docker is None or len(entry["dock"]) >= len(self.docker.targets)
            ):
                continue
            # Pre-filter on Lipinski violations
            if entry["props"]["lipinski_violations"] > self.dock_prefilter_max_violations:
                continue
            to_dock.append(smi)
            to_dock_smiles.append(smi)

        if not to_dock:
            return

        results = self.docker.dock_many(to_dock)
        for smi, result in zip(to_dock_smiles, results):
            self._cache[smi]["dock"] = result  # dict or None (failure)

    # ── Fitness computation ────────────────────────────────────────────

    def _fitness(self, entry: dict[str, Any]) -> tuple[float, dict[str, float], dict[str, Any]]:
        """Compute fitness from cached raw components and current weights."""
        props = entry["props"]
        w = self.weights

        # Normalised components
        parts: dict[str, float] = {}

        # QED: already 0..1
        parts["qed"] = props["qed"]

        # SA: (10 - sa) / 9, so 1.0 (easy) maps to 1.0 and 10.0 (hard) to 0.0
        parts["sa"] = (10.0 - props["sa"]) / 9.0

        # Similarity
        parts["sim"] = entry["sim"]

        # Dock: clip((-E_worst - 5) / 5, 0, 1)
        dock_dict = entry["dock"]
        dock_worst = None
        dock_mean = None
        if dock_dict is not None and len(dock_dict) > 0:
            energies = list(dock_dict.values())
            dock_worst = max(energies)  # weakest (least negative)
            dock_mean = sum(energies) / len(energies)
            parts["dock"] = float(np.clip((-dock_worst - 5.0) / 5.0, 0.0, 1.0))
        else:
            parts["dock"] = 0.0

        # Weighted sum over active components
        active_w = []
        active_c = []
        for name, weight in [("dock", w.dock), ("qed", w.qed), ("sa", w.sa), ("sim", w.sim)]:
            if weight > 0:
                active_w.append(weight)
                active_c.append(parts[name])

        if not active_w:
            fitness = 0.0
        else:
            total_w = sum(active_w)
            fitness = sum(wi * ci for wi, ci in zip(active_w, active_c)) / total_w

        # Penalties
        if props["pains"]:
            fitness *= 0.5

        violations = props["lipinski_violations"]
        if violations > 1:
            fitness *= 0.7 ** (violations - 1)

        # Build raw dict for output
        raw = dict(props)
        raw["dock"] = dock_dict if dock_dict else {}
        raw["dock_worst"] = dock_worst
        raw["dock_mean"] = dock_mean

        return fitness, parts, raw

    # ── Public API ─────────────────────────────────────────────────────

    def score_one(self, smi: str) -> Scored | None:
        """Score a single SMILES. Returns None if the SMILES is invalid."""
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        can = Chem.MolToSmiles(mol)
        self._ensure_cached(mol, can)
        fitness, parts, raw = self._fitness(self._cache[can])
        return Scored(smiles=can, fitness=fitness, parts=parts, raw=raw)

    def score_batch(
        self, smiles_list: list[str], dock: bool = True,
    ) -> list[Scored]:
        """
        Score a batch of SMILES. Triggers docking if enabled.

        Returns list of Scored (only valid molecules included).
        """
        mols_and_smiles: list[tuple[Chem.Mol, str]] = []
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            can = Chem.MolToSmiles(mol)
            self._ensure_cached(mol, can)
            mols_and_smiles.append((mol, can))

        # Batch dock if needed
        if dock:
            self._dock_batch(mols_and_smiles)

        # Compute fitness
        results: list[Scored] = []
        for _mol, can in mols_and_smiles:
            entry = self._cache[can]
            fitness, parts, raw = self._fitness(entry)
            results.append(Scored(smiles=can, fitness=fitness, parts=parts, raw=raw))

        return results

    @property
    def cache_size(self) -> int:
        """Number of molecules in the cache."""
        return len(self._cache)
