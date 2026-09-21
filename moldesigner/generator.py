"""
Genetic molecule generator: BRICS crossover + SMARTS mutations.

All operators are chemistry-aware — every child is valid by construction.
The LLM never generates molecules; this module is the sole source.
"""

from __future__ import annotations

import random
from collections.abc import Generator
from itertools import islice
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem, BRICS, Descriptors, rdChemReactions

from moldesigner.scoring import Scorer, Scored, Weights

# ── Validity filter ────────────────────────────────────────────────────

_REJECT_SMARTS = [
    Chem.MolFromSmarts("[O,S][O,S]"),
    Chem.MolFromSmarts("[#7]~[#7]~[#7]"),
    Chem.MolFromSmarts("[Cl,Br,I][N,O,S]"),
    Chem.MolFromSmarts("C=[C;!R]=C"),
]


def is_valid(mol: Chem.Mol | None) -> bool:
    """
    Check molecule validity: sanitises, single fragment, heavy atoms 12–45,
    MW 150–600, no dangerous substructures.
    """
    if mol is None:
        return False
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return False

    # Single fragment
    frags = Chem.GetMolFrags(mol)
    if len(frags) != 1:
        return False

    # Heavy atom count
    n_heavy = mol.GetNumHeavyAtoms()
    if n_heavy < 12 or n_heavy > 45:
        return False

    # MW
    mw = Descriptors.ExactMolWt(mol)
    if mw < 150 or mw > 600:
        return False

    # Reject dangerous substructures
    for pat in _REJECT_SMARTS:
        if pat is not None and mol.HasSubstructMatch(pat):
            return False

    return True


# ── Mutation SMARTS (all 26 from spec, verified) ──────────────────────

_MUTATION_SMARTS: list[tuple[str, str]] = [
    # Aromatic H substitutions
    ("[cH:1]", "[c:1]F"),
    ("[cH:1]", "[c:1]Cl"),
    ("[cH:1]", "[c:1]C"),
    ("[cH:1]", "[c:1]O"),
    ("[cH:1]", "[c:1]N"),
    ("[cH:1]", "[c:1]OC"),
    ("[cH:1]", "[c:1]C(F)(F)F"),
    ("[cH:1]", "[c:1]C#N"),
    ("[cH:1]", "[c:1]C(=O)N"),
    ("[cH:1]", "[c:1]S(=O)(=O)N"),
    ("[cH:1]", "[c:1]N1CCOCC1"),
    ("[cH:1]", "[c:1]N1CCNCC1"),
    # Aromatic N swaps
    ("[cH:1]", "[n:1]"),
    ("[nH0;R:1]", "[c:1]"),
    # Aliphatic mutations
    ("[CH2:1]", "[O:1]"),
    ("[CH2:1]", "[NH:1]"),
    ("[CH3:1]", "[CH2:1]C"),
    ("[CH3:1]", "[CH2:1]O"),
    ("[OH:1]", "[O:1]C"),
    ("[NH2:1]", "[NH:1]C"),
    ("[NH2:1]", "[NH:1]C(=O)C"),
    ("[NH:1]", "[N:1]C"),
    # Halogen swaps
    ("[F:1]", "[Cl:1]"),
    ("[Cl:1]", "[F:1]"),
    # Functional group swaps
    ("[OH:1]", "[NH2:1]"),
    # Deletion (stripping substituent)
    ("[c:1]-[F,Cl,Br,$([CH3]),$([OH]),$([NH2]),$(C#N)]", "[cH:1]"),
    
    # --- Macro Mutations (Large fragment additions/bioisosteres) ---
    # Swap phenyl with pyridine
    ("c1ccccc1-[*:1]", "n1ccccc1-[*:1]"),
    ("c1ccc([*:1])cc1", "n1ccc([*:1])cc1"),
    # Morpholine addition to aromatic halogen
    ("[c:1]-[F,Cl,Br]", "[c:1]-N1CCOCC1"),
    # Piperazine addition
    ("[c:1]-[F,Cl,Br]", "[c:1]-N1CCNCC1"),
    # Add trifluoromethoxy
    ("[c:1]-[OH]", "[c:1]-OC(F)(F)F"),
    # Direct sulfonyl amide addition
    ("[c:1]-[NH2]", "[c:1]-N(C)S(=O)(=O)C"),
]

# Pre-compile reactions
_MUTATIONS: list[rdChemReactions.ChemicalReaction] = []
for reactant, product in _MUTATION_SMARTS:
    rxn_smarts = f"{reactant}>>{product}"
    rxn = rdChemReactions.ReactionFromSmarts(rxn_smarts)
    if rxn is not None:
        _MUTATIONS.append(rxn)


def mutate(mol: Chem.Mol, rng: random.Random) -> Chem.Mol | None:
    """
    Apply a single random mutation. Returns a valid product or None.

    Shuffles operators, for each tries up to 5 products, keeps the first
    valid one. Converts RunReactants tuple to list before shuffling.
    """
    indices = list(range(len(_MUTATIONS)))
    rng.shuffle(indices)

    for idx in indices:
        rxn = _MUTATIONS[idx]
        try:
            products_tuple = rxn.RunReactants((mol,))
        except Exception:
            continue

        # Convert tuple to list before shuffling (spec §14 gotcha #2)
        products_list = list(products_tuple)
        rng.shuffle(products_list)

        for product_set in products_list[:5]:
            prod = product_set[0]
            try:
                Chem.SanitizeMol(prod)
            except Exception:
                continue
            if is_valid(prod):
                return prod

    return None


_BRICS_CACHE: dict[str, set[str]] = {}

def _decompose_cached(mol: Chem.Mol) -> set[str]:
    smi = Chem.MolToSmiles(mol)
    if smi in _BRICS_CACHE:
        return _BRICS_CACHE[smi]
    try:
        frags = BRICS.BRICSDecompose(mol)
    except Exception:
        frags = set()
    if len(_BRICS_CACHE) > 500:
        _BRICS_CACHE.clear()
    _BRICS_CACHE[smi] = frags
    return frags

def crossover(
    parent1: Chem.Mol, parent2: Chem.Mol, rng: random.Random,
) -> Chem.Mol | None:
    """
    BRICS crossover. Decomposes both parents, unions fragments, builds
    new molecules.

    BRICSBuild has no seedNum argument (spec §14 gotcha #1);
    seed Python global random from our RNG for reproducibility.
    """
    try:
        frags1 = _decompose_cached(parent1)
        frags2 = _decompose_cached(parent2)
    except Exception:
        return None

    all_frags = list(frags1 | frags2)
    if len(all_frags) < 2:
        return None

    rng.shuffle(all_frags)
    subset = all_frags[:6]

    # Parse fragment SMILES to mols
    frag_mols = []
    for f in subset:
        m = Chem.MolFromSmiles(f)
        if m is not None:
            frag_mols.append(m)

    if len(frag_mols) < 2:
        return None

    # Seed Python's global random for BRICSBuild reproducibility
    random.seed(rng.randint(0, 2**31))

    try:
        builder = BRICS.BRICSBuild(frag_mols, onlyCompleteMols=True, maxDepth=2)
        candidates = list(islice(builder, 8))
    except Exception:
        return None

    rng.shuffle(candidates)
    for cand in candidates:
        try:
            Chem.SanitizeMol(cand)
        except Exception:
            continue
        if is_valid(cand):
            return cand

    return None


# ── Morgan fingerprint for diversity ──────────────────────────────────
from rdkit.Chem import rdFingerprintGenerator
_mfp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def _morgan_fp(mol: Chem.Mol):
    return _mfp_gen.GetFingerprint(mol)


def _tanimoto(fp1, fp2) -> float:
    return DataStructs.TanimotoSimilarity(fp1, fp2)


# ── Evolution engine ──────────────────────────────────────────────────

def evolve(
    scorer: Scorer,
    seeds: list[str],
    pop_size: int = 60,
    generations: int = 30,
    crossover_frac: float = 0.5,
    offspring_per_gen: int | None = None,
    diversity_max_sim: float = 0.75,
    seed: int = 0,
) -> Generator[dict[str, Any], None, None]:
    """
    Lazy genetic algorithm generator. Yields one event per generation.

    Callers may change scorer.weights between next() calls;
    the whole population is rescored at the start of each generation.
    """
    if offspring_per_gen is None:
        offspring_per_gen = pop_size

    rng = random.Random(seed)
    np_rng = np.random.RandomState(seed)
    seen: set[str] = set()
    history: list[dict[str, float]] = []

    # ── Initialise population ──────────────────────────────────────
    pop_mols: list[Chem.Mol] = []
    pop_smiles: list[str] = []

    for smi in seeds:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        can = Chem.MolToSmiles(mol)
        if can in seen:
            continue
        if is_valid(mol):
            pop_mols.append(mol)
            pop_smiles.append(can)
            seen.add(can)

    # Fill up with mutants of seeds
    attempts = 0
    max_attempts = pop_size * 20
    while len(pop_mols) < pop_size and attempts < max_attempts:
        attempts += 1
        parent = pop_mols[rng.randint(0, len(pop_mols) - 1)]
        child = mutate(parent, rng)
        if child is None:
            continue
        can = Chem.MolToSmiles(child)
        if can in seen:
            continue
        pop_mols.append(child)
        pop_smiles.append(can)
        seen.add(can)

    # Score initial population
    scored = scorer.score_batch(pop_smiles)
    scored.sort(key=lambda s: s.fitness, reverse=True)

    # Yield generation 0
    gen0_event = _make_event(0, scored, scorer.weights, history)
    history.append({"generation": 0, "best": gen0_event["best"]["fitness"], "mean": gen0_event["mean_fitness"]})
    yield gen0_event

    # ── Main loop ──────────────────────────────────────────────────
    for gen in range(1, generations + 1):
        # Rescore whole population (cheap, uses cache) — applies any weight changes
        pop_smiles_current = [s.smiles for s in scored]
        scored = scorer.score_batch(pop_smiles_current)
        scored.sort(key=lambda s: s.fitness, reverse=True)

        # Build mol lookup
        mol_lookup: dict[str, Chem.Mol] = {}
        for smi in pop_smiles_current:
            m = Chem.MolFromSmiles(smi)
            if m is not None:
                mol_lookup[Chem.MolToSmiles(m)] = m

        # Generate offspring
        offspring_smiles: list[str] = []
        new_count = 0
        off_attempts = 0
        max_off_attempts = offspring_per_gen * 5
        cx_fails = 0  # consecutive crossover failures

        while new_count < offspring_per_gen and off_attempts < max_off_attempts:
            off_attempts += 1

            # Tournament selection (k=3)
            p1 = _tournament(scored, k=3, rng=rng)
            p2 = _tournament(scored, k=3, rng=rng)

            m1 = mol_lookup.get(p1.smiles)
            m2 = mol_lookup.get(p2.smiles)
            if m1 is None or m2 is None:
                continue

            # Crossover or mutation — fall back to mutation if crossover
            # keeps failing (>5 consecutive failures)
            use_cx = rng.random() < crossover_frac and cx_fails < 5
            if use_cx:
                child = crossover(m1, m2, rng)
                if child is None:
                    cx_fails += 1
                    # Fallback: try mutation instead of wasting this attempt
                    child = mutate(m1, rng)
                else:
                    cx_fails = 0
            else:
                child = mutate(m1, rng)
                # 30% chance of second mutation
                if child is not None and rng.random() < 0.3:
                    child2 = mutate(child, rng)
                    if child2 is not None:
                        child = child2

            if child is None:
                continue

            can = Chem.MolToSmiles(child)
            if can in seen:
                continue

            seen.add(can)
            offspring_smiles.append(can)
            new_count += 1

        # Score offspring
        offspring_scored = scorer.score_batch(offspring_smiles)

        # Merge and select survivors with diversity filter
        all_scored = scored + offspring_scored
        all_scored.sort(key=lambda s: s.fitness, reverse=True)

        survivors = _diversity_select(all_scored, pop_size, diversity_max_sim)

        scored = survivors

        # Yield event
        event = _make_event(gen, scored, scorer.weights, history, new_molecules=new_count)
        history.append({"generation": gen, "best": event["best"]["fitness"], "mean": event["mean_fitness"]})
        yield event


def _tournament(
    population: list[Scored], k: int, rng: random.Random,
) -> Scored:
    """Tournament selection: pick k random, return the fittest."""
    contestants = rng.sample(population, min(k, len(population)))
    return max(contestants, key=lambda s: s.fitness)


def _diversity_select(
    sorted_candidates: list[Scored],
    n: int,
    max_sim: float,
) -> list[Scored]:
    """
    Greedy best-first diversity selection.

    Keeps molecules with max Tanimoto ≤ max_sim to those already kept.
    Backfills from leftovers if short.
    """
    kept: list[Scored] = []
    kept_fps: list = []
    leftovers: list[Scored] = []

    for s in sorted_candidates:
        if len(kept) >= n:
            break

        mol = Chem.MolFromSmiles(s.smiles)
        if mol is None:
            continue

        fp = _morgan_fp(mol)

        # Check diversity against kept
        too_similar = False
        for kfp in kept_fps:
            if _tanimoto(fp, kfp) > max_sim:
                too_similar = True
                break

        if too_similar:
            leftovers.append(s)
        else:
            kept.append(s)
            kept_fps.append(fp)

    # Backfill from leftovers if short
    for s in leftovers:
        if len(kept) >= n:
            break
        kept.append(s)

    return kept


def _make_event(
    generation: int,
    scored_pop: list[Scored],
    weights: Weights,
    history: list[dict],
    new_molecules: int = 0,
) -> dict[str, Any]:
    """Build the event dict for one generation."""
    if not scored_pop:
        return {
            "generation": generation,
            "best": None,
            "top": [],
            "mean_fitness": 0.0,
            "new_molecules": new_molecules,
            "weights": weights.as_dict(),
            "history": list(history),
        }

    best = scored_pop[0]
    top = scored_pop[:10]
    mean_f = sum(s.fitness for s in scored_pop) / len(scored_pop)

    def _scored_dict(s: Scored) -> dict:
        return {
            "smiles": s.smiles,
            "fitness": round(s.fitness, 4),
            "parts": {k: round(v, 4) for k, v in s.parts.items()},
            "raw": s.raw,
        }

    return {
        "generation": generation,
        "best": _scored_dict(best),
        "top": [_scored_dict(s) for s in top],
        "mean_fitness": round(mean_f, 4),
        "new_molecules": new_molecules,
        "weights": weights.as_dict(),
        "history": list(history),
    }


def evolve_islands(
    scorer: Scorer,
    seeds: list[str],
    islands: int = 3,
    pop_size_per_island: int = 30,
    generations: int = 30,
    migration_interval: int = 5,
    migration_frac: float = 0.1,
    crossover_frac: float = 0.5,
    diversity_max_sim: float = 0.75,
    seed: int = 0
) -> Generator[dict[str, Any], None, None]:
    """
    Island-model genetic algorithm. Maintains N isolated sub-populations 
    that evolve independently and only communicate during scheduled migration events.
    Yields aggregate events reflecting the global meta-population.
    """
    rng = random.Random(seed)
    
    # Initialize separate island generators
    island_gens = []
    
    # Distribute diversity by shifting seeds slightly per island, if possible
    # We'll just run normal evolve with different RNG seeds per island
    
    island_populations: list[list[Scored]] = [[] for _ in range(islands)]
    
    # We will step standard evolution per island manually
    # Reuse the logic of initialization
    seen: set[str] = set()
    history: list[dict[str, float]] = []

    # Initialize all islands
    for i in range(islands):
        pop_mols = []
        pop_smiles = []
        island_rng = random.Random(seed + i * 100)
        
        for smi in seeds:
            mol = Chem.MolFromSmiles(smi)
            if mol is None: continue
            can = Chem.MolToSmiles(mol)
            if is_valid(mol) and can not in seen:
                pop_mols.append(mol)
                pop_smiles.append(can)
                seen.add(can)
                
        # Fill
        attempts = 0
        while len(pop_mols) < pop_size_per_island and attempts < pop_size_per_island * 20:
            attempts += 1
            if len(pop_mols) == 0: break
            parent = pop_mols[island_rng.randint(0, len(pop_mols) - 1)]
            child = mutate(parent, island_rng)
            if child is None: continue
            can = Chem.MolToSmiles(child)
            if can not in seen:
                pop_mols.append(child)
                pop_smiles.append(can)
                seen.add(can)
                
        scored = scorer.score_batch(pop_smiles)
        scored.sort(key=lambda s: s.fitness, reverse=True)
        island_populations[i] = scored

    # Yield Generation 0 (merge islands for global view)
    all_scored = []
    for pop in island_populations:
        all_scored.extend(pop)
    all_scored.sort(key=lambda s: s.fitness, reverse=True)
    
    # Deduplicate for global event
    global_pop = _diversity_select(all_scored, pop_size_per_island * islands, 0.99)
    event = _make_event(0, global_pop, scorer.weights, history)
    history.append({"generation": 0, "best": event["best"]["fitness"] if event["best"] else 0.0, "mean": event["mean_fitness"]})
    yield event
    
    # Main generational loop across islands
    for gen in range(1, generations + 1):
        total_new_count = 0
        
        for i in range(islands):
            island_rng = random.Random(seed + i * 100 + gen)
            
            # Rescore current island if weights changed
            pop_smiles = [s.smiles for s in island_populations[i]]
            scored = scorer.score_batch(pop_smiles)
            scored.sort(key=lambda s: s.fitness, reverse=True)
            
            mol_lookup: dict[str, Chem.Mol] = {}
            for smi in pop_smiles:
                m = Chem.MolFromSmiles(smi)
                if m is not None:
                    mol_lookup[Chem.MolToSmiles(m)] = m

            offspring_smiles = []
            new_count = 0
            attempts = 0
            cx_fails = 0
            
            while new_count < pop_size_per_island and attempts < pop_size_per_island * 3:
                attempts += 1
                
                # Tournament k=3
                if len(scored) == 0: break
                p1 = _tournament(scored, 3, island_rng)
                p2 = _tournament(scored, 3, island_rng)
                
                m1 = mol_lookup.get(p1.smiles)
                m2 = mol_lookup.get(p2.smiles)
                if m1 is None or m2 is None: continue
                
                if island_rng.random() < crossover_frac and cx_fails < 3:
                    child = crossover(m1, m2, island_rng)
                    if child is None:
                        cx_fails += 1
                        child = mutate(m1, island_rng)
                    else:
                        cx_fails = 0
                else:
                    child = mutate(m1, island_rng)
                    
                if child is None: continue
                can = Chem.MolToSmiles(child)
                if can in seen: continue
                seen.add(can)
                offspring_smiles.append(can)
                new_count += 1
                
            total_new_count += new_count
            offspring_scored = scorer.score_batch(offspring_smiles)
            
            combined = scored + offspring_scored
            combined.sort(key=lambda s: s.fitness, reverse=True)
            island_populations[i] = _diversity_select(combined, pop_size_per_island, diversity_max_sim)
            
        # Migration Event
        if gen % migration_interval == 0 and islands > 1:
            migrate_count = max(1, int(pop_size_per_island * migration_frac))
            migrants: list[list[Scored]] = []
            
            # Extract top N from each island
            for i in range(islands):
                migrants.append(island_populations[i][:migrate_count])
            
            # Ring migration (island i sends to island i+1)
            for i in range(islands):
                received = migrants[i - 1]
                # Replace worst in island with immigrants
                island_populations[i] = island_populations[i][:-migrate_count] + received
                island_populations[i].sort(key=lambda s: s.fitness, reverse=True)

        # Yield global view
        all_scored = []
        for pop in island_populations:
            all_scored.extend(pop)
        all_scored.sort(key=lambda s: s.fitness, reverse=True)
        # Deduplicate
        seen_yield = set()
        global_pop = []
        for s in all_scored:
            if s.smiles not in seen_yield:
                seen_yield.add(s.smiles)
                global_pop.append(s)
                if len(global_pop) >= pop_size_per_island * islands: break
                
        event = _make_event(gen, global_pop, scorer.weights, history, new_molecules=total_new_count)
        history.append({"generation": gen, "best": event["best"]["fitness"] if event["best"] else 0.0, "mean": event["mean_fitness"]})
        yield event

