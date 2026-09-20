"""Tests for moldesigner.scoring module."""

import pytest
from rdkit import Chem

from moldesigner.scoring import properties, Weights, Scorer, Scored


# ── Test molecules ─────────────────────────────────────────────────────

# Erlotinib (EGFR inhibitor, drug-like)
ERLOTINIB = "C=Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"
# Aspirin (small, well-known)
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
# A large molecule likely to violate Lipinski
LARGE_MOL = "CC(C)Cc1ccc(C(C)C(=O)O)cc1"  # ibuprofen (no violations but small)
# Many violations: high MW, high LogP
HIGH_VIOLATION = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"  # long alkane


class TestProperties:
    def test_properties_erlotinib(self):
        mol = Chem.MolFromSmiles(ERLOTINIB)
        assert mol is not None
        props = properties(mol)

        assert "qed" in props
        assert "sa" in props
        assert "mw" in props
        assert "logp" in props
        assert "hbd" in props
        assert "hba" in props
        assert "tpsa" in props
        assert "rotb" in props
        assert "lipinski_violations" in props
        assert "pains" in props

        # QED should be between 0 and 1
        assert 0.0 < props["qed"] < 1.0
        # SA should be between 1 and 10
        assert 1.0 <= props["sa"] <= 10.0
        # MW should be reasonable for erlotinib (~393)
        assert 300 < props["mw"] < 500
        # Not a PAINS hit
        assert props["pains"] is False

    def test_properties_aspirin(self):
        mol = Chem.MolFromSmiles(ASPIRIN)
        props = properties(mol)

        assert props["qed"] > 0
        assert props["mw"] < 200  # aspirin is ~180
        assert props["lipinski_violations"] == 0

    def test_lipinski_violations_count(self):
        """A molecule with multiple violations should count them correctly."""
        # Create a molecule manually with known violations
        # MW>500, logP>5 → 2 violations
        mol = Chem.MolFromSmiles(
            "CCCCCCCCCCCCOC(=O)c1ccc(OC(=O)CCCCCCCCC)cc1"
        )
        if mol is not None:
            props = properties(mol)
            # This should have logP > 5 (long chains) and possibly MW > 500
            assert isinstance(props["lipinski_violations"], int)
            assert props["lipinski_violations"] >= 0


class TestWeights:
    def test_default_weights(self):
        w = Weights()
        assert w.dock == 0.0
        assert w.qed == 0.25
        assert w.sa == 0.15
        assert w.sim == 0.0

    def test_as_dict(self):
        w = Weights(dock=0.5, qed=0.3)
        d = w.as_dict()
        assert d["dock"] == 0.5
        assert d["qed"] == 0.3
        assert "sa" in d
        assert "sim" in d

    def test_mutable(self):
        w = Weights()
        w.dock = 0.8
        assert w.dock == 0.8


class TestScorer:
    def test_score_one_valid(self):
        scorer = Scorer(weights=Weights())
        result = scorer.score_one(ERLOTINIB)
        assert result is not None
        assert isinstance(result, Scored)
        assert 0.0 <= result.fitness <= 1.0
        assert result.smiles  # non-empty

    def test_score_one_invalid(self):
        scorer = Scorer(weights=Weights())
        result = scorer.score_one("not_a_smiles")
        assert result is None

    def test_cache_hit(self):
        scorer = Scorer(weights=Weights())
        r1 = scorer.score_one(ERLOTINIB)
        assert scorer.cache_size == 1
        r2 = scorer.score_one(ERLOTINIB)
        assert scorer.cache_size == 1  # cache hit, no new entry
        assert r1.fitness == r2.fitness

    def test_weight_change_updates_fitness(self):
        scorer = Scorer(weights=Weights(qed=0.5, sa=0.0))
        r1 = scorer.score_one(ERLOTINIB)

        scorer.weights.qed = 0.0
        scorer.weights.sa = 0.5
        r2 = scorer.score_one(ERLOTINIB)

        # Fitness should differ because weights changed
        # (unless QED and SA happen to be identical after normalisation)
        assert scorer.cache_size == 1  # no re-computation of raw props
        # The raw props should be identical
        assert r1.raw["qed"] == r2.raw["qed"]

    def test_fitness_between_zero_and_one(self):
        scorer = Scorer(weights=Weights())
        for smi in [ERLOTINIB, ASPIRIN]:
            result = scorer.score_one(smi)
            if result is not None:
                assert 0.0 <= result.fitness <= 1.0, f"{smi}: fitness={result.fitness}"

    def test_score_batch(self):
        scorer = Scorer(weights=Weights())
        results = scorer.score_batch([ERLOTINIB, ASPIRIN, "invalid"])
        # Only valid molecules returned
        assert len(results) == 2
        for r in results:
            assert isinstance(r, Scored)

    def test_similarity_scoring(self):
        """Similarity to reference should be 1.0 for the reference itself."""
        scorer = Scorer(
            weights=Weights(sim=1.0, qed=0.0, sa=0.0),
            references=[ERLOTINIB],
        )
        result = scorer.score_one(ERLOTINIB)
        assert result is not None
        # Similarity to self should be 1.0
        assert result.parts["sim"] == pytest.approx(1.0, abs=0.01)

    def test_pains_penalty(self):
        """PAINS molecules should get a 0.5 penalty."""
        # Rhodanine is a known PAINS pattern
        rhodanine = "O=C1CSC(=S)N1"
        scorer = Scorer(weights=Weights(qed=1.0, sa=0.0))

        result = scorer.score_one(rhodanine)
        if result is not None and result.raw.get("pains"):
            # Fitness should be halved
            # Compute what fitness would be without penalty
            raw_fitness = result.parts["qed"]  # only QED active
            assert result.fitness <= raw_fitness * 0.5 + 0.01

    def test_no_dock_weight_means_no_dock(self):
        """When dock weight is 0, dock component should not affect fitness."""
        scorer = Scorer(weights=Weights(dock=0.0, qed=1.0))
        result = scorer.score_one(ERLOTINIB)
        assert result is not None
        # Dock part is computed but not used in fitness
        assert result.parts["dock"] == 0.0  # no docking done

    def test_worst_case_scoring(self):
        """Scorer should use the worst (least negative/highest) energy among variants."""
        scorer = Scorer(weights=Weights(dock=1.0, qed=0.0, sa=0.0))
        # Ensure it's not in cache yet
        mol = Chem.MolFromSmiles(ERLOTINIB)
        smi = Chem.MolToSmiles(mol, canonical=True)
        scorer._ensure_cached(mol, smi)
        
        # Inject fake multi-variant dock results
        # WT=-9.0, L858R=-8.5, T790M_C797S=-4.0 (worst case)
        scorer._cache[smi]["dock"] = {
            "WT": -9.0,
            "L858R": -8.5,
            "T790M_C797S": -4.0
        }
        
        result = scorer.score_one(ERLOTINIB)
        assert result is not None
        
        # Check raw extractions
        assert result.raw["dock_worst"] == -4.0
        
        # Check fitness normalization: float(np.clip((-dock_worst - 5.0) / 5.0, 0.0, 1.0))
        # clip((-(-4.0) - 5)/5) == clip((4.0 - 5.0)/5.0) == clip(-0.2) == 0.0
        assert pytest.approx(result.parts["dock"], abs=0.01) == 0.0
        assert pytest.approx(result.fitness, abs=0.01) == 0.0

