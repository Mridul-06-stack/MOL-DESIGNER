# MolDesigner: Adversarial AI Drug Discovery 🧬

MolDesigner is an advanced, closed-loop Artificial Intelligence drug discovery platform built to programmatically evolve molecular compounds capable of circumventing biological tumor resistance (such as the lung cancer EGFR mutations).

By combining an island-model continuous genetic algorithm (GA) with an automated "red-team" protein mutation simulator, MolDesigner acts as a continuous self-playing game engine against cancer resistance.

---

## 🛠 Architecture Flow & Block Diagram

```mermaid
graph TD
    %% Define Styles
    classDef react fill:#00d8ff,stroke:#000,stroke-width:2px,color:#000;
    classDef fastAPI fill:#009688,stroke:#000,stroke-width:2px,color:#fff;
    classDef ai fill:#673ab7,stroke:#000,stroke-width:2px,color:#fff;
    classDef bio fill:#e91e63,stroke:#000,stroke-width:2px,color:#fff;

    subgraph Frontend [React+Vite UI Application]
        UI[Glassmorphic React Dashboard]:::react
        Viewer[3Dmol.js Molecular Visualizer]:::react
    end

    subgraph Backend [Python FastAPI Server]
        API[POST /evolve (Server-Sent Events)]:::fastAPI
        MolAPI[GET /molblock (3D Coordinates)]:::fastAPI
    end

    subgraph Core [MolDesigner Intelligent Core]
        Gen[Island-Model Genetic Generator]:::ai
        Mutate((Structural Mutator\n& Crossover)):::ai
        
        Dock[Docker Interface]:::bio
        Vina[AutoDock Vina / Meeko]:::bio
        RDKit[RDKit Fallback Scoring]:::bio
        
        Scanner[Red-Team Escape Scanner]:::bio
    end

    %% Flow Dynamics
    UI -- "Config (Islands, Target Panel)" --> API
    UI -- "Fetch 3D Mapping" --> MolAPI
    
    API -- "Initiate & Stream Gen Events" --> Gen
    
    Gen -- "Spawn Chemical Variants" --> Mutate
    Mutate -- "Produces SMILES Array" --> Dock
    
    Dock -- "Evaluates Binding Energies" --> Vina
    Dock -- "Fallback Fast Evaluation" --> RDKit
    
    Vina -. "Returns Binding Severities" .-> Gen
    RDKit -. "Returns Multi-variant Scores" .-> Gen
    
    %% Adversarial Loop
    Gen -- "Sends Best Surviving Candidate" --> Scanner
    Scanner -- "Identifies Weakness & Injects New Protein Mutation" --> Dock
    
    Gen -- "Yields Live SSE Stream" --> UI
```

---

## 🚀 The AI adversarial Loop

MolDesigner implements a Continuous Evolution paradigm:
1. **The Generator (White Team):** Rapidly iterates molecular geometries via fragment-swapping and bioisostere crossover, identifying a molecule capable of neutralizing modern protein mutations (e.g. wild-type and L858R combined).
2. **The Docker (Physics Engine):** Ensures all surviving molecules adhere to multi-variant structural constraints without falling victim to localized optimums.
3. **The Simulator (Red Team):** Ingests the White Team's best surviving drug candidate, mathematically identifying a resistance loophole, and permanently injecting a new synthetic protein-mutant target into the Docker configuration to forcibly disrupt the drug's efficacy and restart the cycle.

---

## 💻 Requirements & Installation

MolDesigner fundamentally breaks down into a heavy numerical backend and a visual frontend layer.

### 1. Python Backend Dependencies

You must install standard high-capacity cheminformatics processing.
**Tooling Requirements:** Python 3.10+, pip

```bash
# Core Chemistry AI dependencies
pip install rdkit-pypi scipy numpy pytest

# Standalone Server requirements
pip install fastapi uvicorn pydantic

# (Optional: Only if employing physical binding instead of surrogate models)
pip install vina meeko
```

### 2. React Frontend Dependencies

The UX drives real-time analysis through Server-Sent Event streaming.
**Tooling Requirements:** Node.js (v18+), npm

```bash
cd ui
npm install
```

---

## ⚡ Execution Commands

MolDesigner is designed with full multi-stage capabilities. Start both modules concurrently to use the unified console.

#### 1. Boot up the API Engine
```bash
uvicorn server.api:app --host 0.0.0.0 --port 8000
```
This initializes the FastAPI sub-layer waiting for commands from the dashboard.

#### 2. Launch the Control UI
```bash
cd ui
npm run dev
```
Navigate to `http://localhost:5173` to access the visual real-time generator and 3D molecular renderer!
