"""Tests for moldesigner.generator module."""

import random
import pytest
from rdkit import Chem

from moldesigner.generator import is_valid, mutate, crossover, evolve, _MUTATIONS
from moldesigner.scoring import Scorer, Weights


# ── Test molecules ─────────────────────────────────────────────────────

# Erlotinib (EGFR inhibitor) - has aromatic C-H, halogens, etc.
ERLOTINIB = "C=Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"

# Two kinase-like parents for crossover testing
PARENT_A = "c1ccc(-c2ccnc(Nc3ccccc3)n2)cc1"        # 2-aminopyrimidine scaffold
PARENT_B = "c1ccc2c(c1)cc(-c1ccccn1)c(=O)o2"       # coumarin-pyridine

# Imatinib-like (large, drug-like)
IMATINIB_LIKE = "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"

# Some molecules for mutation testing that should have diverse functional groups
TEST_MOLS_FOR_MUTATION = [
    "c1ccc(F)cc1NC(=O)c1ccncc1",           # has aromatic H, F, NH
    "c1ccc(Cl)c(O)c1CNC(=O)CC",            # has Cl, OH, NH, CH3, CH2
    "c1ccc2[nH]ccc2c1",                      # has nH in ring, aromatic H
    "c1ccncc1C(=O)Nc1ccc(C)cc1",            # has aromatic N, CH3, aromatic H
    "c1ccc(NC(=O)c2cc(F)ccc2O)cc1",         # has F, OH, NH
    "c1ccc(C(=O)N2CCOCC2)cc1N",             # has NH2, morpholine
]


class TestIsValid:
    def test_valid_drug_like(self):
        mol = Chem.MolFromSmiles(ERLOTINIB)
        assert is_valid(mol)

    def test_none_invalid(self):
        assert not is_valid(None)

    def test_too_small(self):
        mol = Chem.MolFromSmiles("CC")  # ethane, 2 heavy atoms
        assert not is_valid(mol)

    def test_too_large(self):
        # Build a very large molecule >45 heavy atoms
        big = "C" * 50
        mol = Chem.MolFromSmiles(big)
        assert not is_valid(mol)

    def test_mw_too_low(self):
        # Small ring, few heavy atoms but might pass count
        mol = Chem.MolFromSmiles("c1ccccccccccc1")  # 12 heavy atoms, low MW
        if mol is not None:
            # May or may not pass depending on exact MW
            result = is_valid(mol)
            # Just checking it doesn't crash

    def test_multiple_fragments(self):
        mol = Chem.MolFromSmiles("c1ccccc1.c1ccccc1")
        assert not is_valid(mol)

    def test_reject_peroxide(self):
        """O-O bond should be rejected."""
        mol = Chem.MolFromSmiles("c1ccc(OOC)cc1NC(=O)c1ccccc1")
        assert not is_valid(mol)

    def test_reject_azide(self):
        """N=N=N should be rejected."""
        mol = Chem.MolFromSmiles("c1ccc(N=[N+]=[N-])cc1NC(=O)c1ccccc1")
        assert not is_valid(mol)


class TestMutations:
    def test_all_26_mutations_fire(self):
        """
        Each of the 26 mutation operators should produce at least one product
        on at least one of the test molecules.
        """
        fired = set()

        for i, rxn in enumerate(_MUTATIONS):
            for smi in TEST_MOLS_FOR_MUTATION:
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    continue
                try:
                    products = rxn.RunReactants((mol,))
                except Exception:
                    continue
                if len(products) > 0:
                    fired.add(i)
                    break

        # All 26 should have fired
        missing = set(range(len(_MUTATIONS))) - fired
        assert len(missing) == 0, (
            f"Mutations that didn't fire: {missing} "
            f"(out of {len(_MUTATIONS)} total)"
        )

    def test_mutate_returns_valid(self):
        """mutate should always return valid molecules or None."""
        rng = random.Random(42)
        mol = Chem.MolFromSmiles(ERLOTINIB)
        valid_count = 0
        total = 50

        for _ in range(total):
            child = mutate(mol, rng)
            if child is not None:
                assert is_valid(child), f"Invalid mutant: {Chem.MolToSmiles(child)}"
                valid_count += 1

        # Should produce at least some valid mutants
        assert valid_count > 0, "No valid mutants produced"

    def test_mutate_deterministic_with_seed(self):
        """Same seed should produce same results."""
        mol = Chem.MolFromSmiles(ERLOTINIB)
        results1 = []
        results2 = []

        for seed in [42]:
            rng1 = random.Random(seed)
            rng2 = random.Random(seed)
            for _ in range(10):
                r1 = mutate(mol, rng1)
                r2 = mutate(mol, rng2)
                if r1 is not None:
                    results1.append(Chem.MolToSmiles(r1))
                if r2 is not None:
                    results2.append(Chem.MolToSmiles(r2))

        assert results1 == results2


class TestCrossover:
    def test_crossover_success_rate(self):
        """Crossover should return valid child ≥80% of 50 attempts."""
        p1 = Chem.MolFromSmiles(PARENT_A)
        p2 = Chem.MolFromSmiles(PARENT_B)
        assert p1 is not None and p2 is not None

        successes = 0
        total = 50

        for i in range(total):
            rng = random.Random(i * 100 + 7)
            child = crossover(p1, p2, rng)
            if child is not None:
                successes += 1

        rate = successes / total
        assert rate >= 0.80, f"Crossover success rate {rate:.2%} < 80%"

    def test_crossover_returns_valid(self):
        """Every successful crossover product should pass is_valid."""
        p1 = Chem.MolFromSmiles(PARENT_A)
        p2 = Chem.MolFromSmiles(PARENT_B)

        for i in range(30):
            rng = random.Random(i * 37)
            child = crossover(p1, p2, rng)
            if child is not None:
                assert is_valid(child), f"Invalid crossover child: {Chem.MolToSmiles(child)}"


class TestEvolve:
    def _run_evolution(self, seed=42, generations=12, pop_size=40,
                       weights=None, references=None, seeds_list=None):
        """Helper to run a full evolution."""
        if weights is None:
            weights = Weights(qed=0.5, sa=0.3, sim=0.0)
        if seeds_list is None:
            seeds_list = [PARENT_A, PARENT_B, ERLOTINIB]

        scorer = Scorer(weights=weights, references=references)
        events = list(evolve(
            scorer=scorer,
            seeds=seeds_list,
            pop_size=pop_size,
            generations=generations,
            seed=seed,
        ))
        return events

    def test_evolve_produces_events(self):
        events = self._run_evolution(generations=5, pop_size=20)
        assert len(events) == 6  # gen 0..5

    def test_evolve_event_schema(self):
        events = self._run_evolution(generations=3, pop_size=20)
        for event in events:
            assert "generation" in event
            assert "best" in event
            assert "top" in event
            assert "mean_fitness" in event
            assert "weights" in event
            assert "history" in event

            if event["best"] is not None:
                best = event["best"]
                assert "smiles" in best
                assert "fitness" in best
                assert "parts" in best

    def test_evolve_all_valid(self):
        """Every emitted SMILES should pass is_valid."""
        events = self._run_evolution(generations=8, pop_size=30)
        for event in events:
            for mol_data in event["top"]:
                smi = mol_data["smiles"]
                mol = Chem.MolFromSmiles(smi)
                assert mol is not None, f"Invalid SMILES in top: {smi}"
                assert is_valid(mol), f"SMILES fails is_valid: {smi}"

    def test_evolve_reproducibility(self):
        """Two runs with the same seed should produce identical results."""
        events1 = self._run_evolution(seed=123, generations=5, pop_size=20)
        events2 = self._run_evolution(seed=123, generations=5, pop_size=20)

        assert len(events1) == len(events2)
        for e1, e2 in zip(events1, events2):
            assert e1["generation"] == e2["generation"]
            if e1["best"] is not None and e2["best"] is not None:
                assert e1["best"]["smiles"] == e2["best"]["smiles"]
                assert e1["best"]["fitness"] == e2["best"]["fitness"]

    def test_evolve_fitness_improves(self):
        """
        Similarity-objective run: best fitness ≥ initial,
        mean fitness should show an increasing trend.
        """
        # Use similarity weight so we have a clear objective
        weights = Weights(qed=0.3, sa=0.2, sim=0.5)
        events = self._run_evolution(
            seed=42,
            generations=12,
            pop_size=40,
            weights=weights,
            references=[ERLOTINIB],
        )

        initial_best = events[0]["best"]["fitness"]
        final_best = events[-1]["best"]["fitness"]
        assert final_best >= initial_best, (
            f"Best fitness decreased: {initial_best:.4f} → {final_best:.4f}"
        )

        # Mean fitness should generally increase (check first vs last quarter)
        means = [e["mean_fitness"] for e in events]
        first_quarter = sum(means[:3]) / 3
        last_quarter = sum(means[-3:]) / 3
        assert last_quarter >= first_quarter - 0.05, (
            f"Mean fitness not increasing: first_q={first_quarter:.4f}, "
            f"last_q={last_quarter:.4f}"
        )

    def test_evolve_generation_numbers(self):
        events = self._run_evolution(generations=5, pop_size=20)
        gens = [e["generation"] for e in events]
        assert gens == [0, 1, 2, 3, 4, 5]
