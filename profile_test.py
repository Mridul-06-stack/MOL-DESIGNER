"""Profile evolve with test_evolve_fitness_improves params."""
import time
from moldesigner.scoring import Scorer, Weights
from moldesigner.generator import evolve

ERLOTINIB = "C=Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"
PARENT_A = "c1ccc(-c2ccnc(Nc3ccccc3)n2)cc1"
PARENT_B = "c1ccc2c(c1)cc(-c1ccccn1)c(=O)o2"

w = Weights(qed=0.3, sa=0.2, sim=0.5)
s = Scorer(weights=w, references=[ERLOTINIB])
seeds = [PARENT_A, PARENT_B, ERLOTINIB]

t0 = time.time()
gen = evolve(scorer=s, seeds=seeds, pop_size=40, generations=12, seed=42)
total_new = 0

for i, ev in enumerate(gen):
    t1 = time.time()
    g = ev["generation"]
    best = ev["best"]["fitness"] if ev["best"] else 0
    new = ev["new_molecules"]
    total_new += new
    print(f"gen {g}: {t1-t0:.1f}s, best={best:.4f}, new={new}, total_seen={len(s._cache)}")
    t0 = t1

print(f"Done! Total new: {total_new}")
