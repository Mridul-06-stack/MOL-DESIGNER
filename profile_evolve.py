"""Quick profile of evolve performance."""
import time
from moldesigner.scoring import Scorer, Weights
from moldesigner.generator import evolve

w = Weights(qed=0.5, sa=0.3, sim=0.0)
s = Scorer(weights=w)
seeds = [
    "c1ccc(-c2ccnc(Nc3ccccc3)n2)cc1",
    "c1ccc2c(c1)cc(-c1ccccn1)c(=O)o2",
    "C=Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1",
]

t0 = time.time()
gen = evolve(scorer=s, seeds=seeds, pop_size=20, generations=3, seed=42)
for ev in gen:
    t1 = time.time()
    g = ev["generation"]
    best = ev["best"]["fitness"]
    new = ev["new_molecules"]
    print(f"gen {g}: {t1-t0:.1f}s, best={best:.4f}, new={new}")
    t0 = t1

print("Done!")
