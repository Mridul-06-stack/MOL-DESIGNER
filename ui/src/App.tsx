import { useState, useEffect, useRef } from 'react'

declare global {
  interface Window {
    $3Dmol?: any
  }
}

interface ScoredBest {
  smiles: string
  fitness: number
  parts: {
    qed?: number
    sa?: number
    sim?: number
    dock?: number
    [key: string]: number | undefined
  }
  raw?: {
    qed?: number
    sa?: number
    mw?: number
    logp?: number
    hbd?: number
    hba?: number
    tpsa?: number
    rotb?: number
    lipinski_violations?: number
    pains?: boolean
    dock?: Record<string, number>
    dock_worst?: number | null
    dock_mean?: number | null
    [key: string]: any
  }
}

interface GenerationEvent {
  generation: number
  best: ScoredBest | null
  top: ScoredBest[]
  mean_fitness: number
  new_molecules: number
  weights: Record<string, number>
  history: Array<{ generation: number; best: number; mean: number }>
}

export default function App() {
  const [evolving, setEvolving] = useState(false)
  const [generations, setGenerations] = useState<GenerationEvent[]>([])
  const [currentBest, setCurrentBest] = useState<ScoredBest | null>(null)
  const [statusText, setStatusText] = useState<string>('')
  const viewerRef = useRef<any>(null)

  const [config, setConfig] = useState({
    generations: 20,
    islands: 3,
    pop_size: 30,
    dock_weight: 0.5,
    seed_smiles: "C=Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1",
    preset: "erlotinib"
  })

  // Start evolution request
  const handleStart = async () => {
    setEvolving(true)
    setGenerations([])
    setCurrentBest(null)
    setStatusText('Connecting to evolution engine...')

    try {
      // Stream generation events via SSE
      const res = await fetch('http://localhost:8000/api/evolve', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ ...config, targets: ["WT", "L858R", "T790M_C797S"] })
      })

      if (!res.ok) {
        setStatusText(`Server error: ${res.statusText}`)
        setEvolving(false)
        return
      }

      if (!res.body) {
        setStatusText('No response stream available')
        setEvolving(false)
        return
      }

      const reader = res.body.getReader()
      const decoder = new TextDecoder()
      let buffer = ''
      setStatusText('Simulating evolution against resistance panel...')

      while (true) {
        const { value, done } = await reader.read()
        if (done) break

        buffer += decoder.decode(value, { stream: true })
        const parts = buffer.split('\n\n')
        buffer = parts.pop() || ''

        for (const part of parts) {
          if (part.startsWith('data: ')) {
            try {
              const data: GenerationEvent = JSON.parse(part.slice(6))
              setGenerations(prev => [...prev, data])
              if (data.best) {
                setCurrentBest(data.best)
                render3D(data.best.smiles)
              }
            } catch (e) {
              console.error('JSON Parse error', e)
            }
          }
        }
      }
      setStatusText('Evolution run complete!')
    } catch (err: any) {
      console.error('Evolution error:', err)
      setStatusText(`Connection failed: Make sure FastAPI backend is running on port 8000.`)
    } finally {
      setEvolving(false)
    }
  }

  const render3D = async (smiles: string) => {
    try {
      const res = await fetch(`http://localhost:8000/api/molblock?smiles=${encodeURIComponent(smiles)}`)
      if (!res.ok) return
      const sdf = await res.text()
      if (!sdf) return

      if (typeof window !== 'undefined' && window.$3Dmol) {
        const container = document.getElementById("viewer3d")
        if (!container) return
        container.innerHTML = ""

        const viewer = window.$3Dmol.createViewer(container, {
          defaultcolors: window.$3Dmol.rasmolElementColors
        })
        viewerRef.current = viewer
        viewer.addModel(sdf, "sdf")
        viewer.setStyle({}, { stick: { radius: 0.15 }, sphere: { scale: 0.25 } })
        viewer.zoomTo()
        viewer.render()
      }
    } catch (e) {
      console.error('3D rendering error:', e)
    }
  }

  // Effect to load 3dmol script
  useEffect(() => {
    if (!document.getElementById("3dmol-script")) {
      const script = document.createElement("script")
      script.id = "3dmol-script"
      script.src = "https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"
      script.async = true
      document.body.appendChild(script)
    }
  }, [])

  const worstEnergy = currentBest?.raw?.dock_worst != null
    ? currentBest.raw.dock_worst.toFixed(2)
    : '--'

  return (
    <div className="dashboard-container">
      <header>
        <div>
          <h1>MolDesigner</h1>
          <p style={{ color: 'var(--text-secondary)', fontSize: '0.875rem', marginTop: '0.25rem' }}>
            Closed-Loop Adversarial AI Drug Discovery & EGFR Resistance Overcoming
          </p>
        </div>
        {evolving && (
          <div className="pill pulse" style={{ background: 'rgba(59,130,246,0.2)', color: 'var(--accent-blue)' }}>
            Running AI Evolution...
          </div>
        )}
      </header>

      <aside className="sidebar glass-panel">
        <h2 style={{ fontSize: '1.2rem', marginBottom: '1.5rem', paddingBottom: '0.5rem', borderBottom: '1px solid rgba(255,255,255,0.1)' }}>
          Simulation Parameters
        </h2>

        <div className="form-group">
          <label>Starting Molecule Preset</label>
          <select
            value={config.preset}
            onChange={e => {
              const val = e.target.value
              let sm = config.seed_smiles
              if (val === 'erlotinib') sm = "C=Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"
              else if (val === 'aspirin') sm = "CC(=O)Oc1ccccc1C(=O)O"
              else if (val === 'caffeine') sm = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
              else if (val === 'penicillin') sm = "CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C"
              else if (val === 'imatinib') sm = "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C"
              else if (val === 'custom') sm = ""

              setConfig({ ...config, preset: val, seed_smiles: sm })
            }}
          >
            <option value="erlotinib">Erlotinib (EGFR - Lung Cancer)</option>
            <option value="imatinib">Imatinib (BCR-ABL - Gleevec)</option>
            <option value="penicillin">Penicillin G (Antibiotic)</option>
            <option value="aspirin">Aspirin (Anti-inflammatory)</option>
            <option value="caffeine">Caffeine (Adenosine Antagonist)</option>
            <option value="custom">Custom SMILES</option>
          </select>
        </div>

        <div className="form-group">
          <label>Starting SMILES Seed</label>
          <input
            type="text"
            value={config.seed_smiles}
            onChange={e => setConfig({ ...config, seed_smiles: e.target.value, preset: 'custom' })}
            style={{ fontFamily: 'monospace', fontSize: '0.75rem' }}
          />
        </div>

        <div className="form-group">
          <label>Generations</label>
          <input
            type="number"
            min={1}
            max={100}
            value={config.generations}
            onChange={e => setConfig({ ...config, generations: Math.max(1, parseInt(e.target.value) || 1) })}
          />
        </div>

        <div className="form-group">
          <label>Island Pools (Migration Model)</label>
          <input
            type="number"
            min={1}
            max={8}
            value={config.islands}
            onChange={e => setConfig({ ...config, islands: Math.max(1, parseInt(e.target.value) || 1) })}
          />
        </div>

        <div className="form-group">
          <label>Population Size</label>
          <input
            type="number"
            min={10}
            max={200}
            value={config.pop_size}
            onChange={e => setConfig({ ...config, pop_size: Math.max(10, parseInt(e.target.value) || 10) })}
          />
        </div>

        <button
          className="btn-primary"
          onClick={handleStart}
          disabled={evolving}
          style={{ marginTop: '1rem', opacity: evolving ? 0.6 : 1 }}
        >
          {evolving ? 'Evolving Chemistry...' : 'Launch Generator'}
        </button>

        {statusText && (
          <p style={{ fontSize: '0.8rem', color: 'var(--text-secondary)', marginTop: '0.5rem', textAlign: 'center' }}>
            {statusText}
          </p>
        )}

        {currentBest && (
          <div style={{ marginTop: '1.5rem', padding: '1rem', background: 'rgba(0,0,0,0.3)', borderRadius: '8px', border: '1px solid rgba(255,255,255,0.08)' }}>
            <h3 style={{ fontSize: '0.85rem', color: 'var(--text-secondary)', marginBottom: '0.5rem', textTransform: 'uppercase', letterSpacing: '0.05em' }}>
              Current Best Scaffold
            </h3>
            <p style={{ wordBreak: 'break-all', fontSize: '0.8rem', fontFamily: 'monospace', color: '#60a5fa' }}>
              {currentBest.smiles}
            </p>
            {currentBest.raw && (
              <div style={{ marginTop: '0.75rem', fontSize: '0.75rem', color: 'var(--text-secondary)', display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '0.25rem' }}>
                <div>MW: {currentBest.raw.mw?.toFixed(1)} Da</div>
                <div>LogP: {currentBest.raw.logp?.toFixed(2)}</div>
                <div>QED: {currentBest.raw.qed?.toFixed(2)}</div>
                <div>SA: {currentBest.raw.sa?.toFixed(2)}</div>
              </div>
            )}
          </div>
        )}
      </aside>

      <main className="main-content">
        <div className="stats-grid">
          <div className="stat-card glass-panel">
            <h3>Generations Completed</h3>
            <div className="value">{generations.length} / {config.generations}</div>
          </div>
          <div className="stat-card glass-panel">
            <h3>Best Composite Fitness</h3>
            <div className="value" style={{ color: currentBest && currentBest.fitness > 0.6 ? 'var(--success-color)' : 'white' }}>
              {currentBest ? (currentBest.fitness * 100).toFixed(1) + '%' : '--'}
            </div>
          </div>
          <div className="stat-card glass-panel">
            <h3>Worst-Case Energy (kcal/mol)</h3>
            <div className="value" style={{ color: currentBest?.raw?.dock_worst && currentBest.raw.dock_worst < -8 ? 'var(--success-color)' : 'white' }}>
              {worstEnergy}
            </div>
          </div>
        </div>

        <div className="glass-panel" style={{ display: 'flex', flexDirection: 'column', padding: '1.5rem', minHeight: '520px' }}>
          <h2 style={{ fontSize: '1.1rem', marginBottom: '1rem', display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
            <span>Interactive 3D Molecular Conformer</span>
            <span style={{ fontSize: '0.8rem', color: 'var(--text-secondary)' }}>Powered by 3Dmol.js & RDKit MMFF</span>
          </h2>
          <div id="viewer3d" className="viewer-container" style={{ position: 'relative', width: '100%', height: '450px' }}>
            {!currentBest && !evolving && (
              <div style={{ color: 'var(--text-secondary)', textAlign: 'center', padding: '2rem' }}>
                Click <strong>"Launch Generator"</strong> on the left panel to begin evolution and render 3D conformers.
              </div>
            )}
            {evolving && !currentBest && (
              <div className="loading-spinner"></div>
            )}
          </div>
        </div>
      </main>
    </div>
  )
}

