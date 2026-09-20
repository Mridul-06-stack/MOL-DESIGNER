# MolDesigner

Resistance-aware molecule design with a red-team escape loop.

## What this is

An interactive design loop that evolves drug-like molecules against a panel of protein
variants (wild-type + resistance mutants), then attacks its own best candidates with
predicted escape mutations and redesigns against the ones that break them.

**Honest novelty claim:** No single component is new. We did not find an open, interactive
tool coupling a molecule generator to an automated escape-mutation scanner in a closed loop.
That is a search result, not proof of absence.

**Not:** a drug-discovery product, a validated predictor, or a source of drug candidates.
Outputs are computational hypotheses.

> **Docking scores are proxies** with known weak correlation to real binding affinity.
> Nothing here is experimentally validated.

## Quick start

```bash
# Install
pip install -e ".[dev]"

# Run headless (Phase 1, no docking)
python cli/run_headless.py \
    --seeds "c1ccc2c(c1)cc(=O)oc2" "c1ccc(-c2ccncc2)cc1" \
    --generations 12 --pop-size 40 --seed 42

# Tests
pytest tests/ -v
```

## Phases

1. **Core** — scoring + generator + CLI (no docking)
2. **Docking + validation** — AutoDock Vina panel
3. **EGFR panel** — multi-variant worst-case scoring
4. **Generator strength** — larger moves, island populations
5. **Escape scanner** — in-silico mutation scanning
6. **Server + UI** — FastAPI + React dashboard
7. **Stretch** — GNINA, AiZynthFinder, agent layer
