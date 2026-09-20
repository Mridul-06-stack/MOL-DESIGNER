import { useState, useEffect, useRef } from 'react'

export default function App() {
  const [evolving, setEvolving] = useState(false)
  const [generations, setGenerations] = useState([])
  const [currentBest, setCurrentBest] = useState(null)
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

    // Trigger backend POST
    await fetch('http://localhost:8000/api/evolve', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ ...config, targets: ["WT", "L858R", "T790M_C797S"] })
    }).catch(e => console.error(e))

    // Because FastAPI doesn't easily support POST SSE from browser EventSource native,
    // we use a fetch reader to parse the SSE stream manually.
    const res = await fetch('http://localhost:8000/api/evolve', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ ...config, targets: ["WT", "L858R", "T790M_C797S"] })
    })

    const reader = res.body.getReader()
    const decoder = new TextDecoder()
    let buffer = ''

    while (true) {
      const { value, done } = await reader.read()
      if (done) break

      buffer += decoder.decode(value, { stream: true })

      const parts = buffer.split('\n\n')
      buffer = parts.pop() || ''

      for (const part of parts) {
        if (part.startsWith('data: ')) {
          try {
            const data = JSON.parse(part.slice(6))
            setGenerations(prev => [...prev, data])
            setCurrentBest(data.best)

            // Render 3D if best exists
            if (data.best && typeof window.$3Dmol !== 'undefined') {
              render3D(data.best.smiles)
            }
          } catch (e) {
            console.error('JSON Parse error', e)
          }
        }
      }
    }
    setEvolving(false)
  }

  const render3D = async (smiles) => {
    try {
      const res = await fetch(`http://localhost:8000/api/molblock?smiles=${encodeURIComponent(smiles)}`)
      const sdf = await res.text()

      const viewer = window.$3Dmol.createViewer("viewer3d", { defaultcolors: window.$3Dmol.rasmolElementColors })
      viewer.addModel(sdf, "sdf")
      viewer.setStyle({}, { stick: { radius: 0.15 }, sphere: { scale: 0.25 } })
      viewer.zoomTo()
      viewer.render()
    } catch (e) {
      console.error(e)
    }
  }

  // Effect to load 3dmol script
  useEffect(() => {
    const script = document.createElement("script")
    script.src = "https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"
    script.async = true
    document.body.appendChild(script)
  }, [])

  return (
    <div className="dashboard-container">
      <header>
        <h1>MolDesigner</h1>
        {evolving && <div className="pill pulse" style={{ background: 'rgba(59,130,246,0.2)', color: 'var(--accent-blue)' }}>Running AI Pipeline...</div>}
      </header>

      <aside className="sidebar glass-panel">
        <h2 style={{ fontSize: '1.2rem', marginBottom: '1.5rem', paddingBottom: '0.5rem', borderBottom: '1px solid rgba(255,255,255,0.1)' }}>Simulation Params</h2>

        <div className="form-group">
          <label>Choose Starting Molecule</label>
          <select
            value={config.preset}
            onChange={e => {
              const val = e.target.value;
              let sm = config.seed_smiles;
              if (val === 'erlotinib') sm = "C=Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"
              else if (val === 'aspirin') sm = "CC(=O)Oc1ccccc1C(=O)O"
              else if (val === 'caffeine') sm = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
              else if (val === 'penicillin') sm = "CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C"
              else if (val === 'imatinib') sm = "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C"
              else if (val === 'custom') sm = ""

              setConfig({ ...config, preset: val, seed_smiles: sm })
            }}
          >
            <option value="erlotinib">Erlotinib (Lung Cancer)</option>
            <option value="imatinib">Imatinib (Leukemia - Gleevec)</option>
            <option value="penicillin">Penicillin G (Antibiotic)</option>
            <option value="aspirin">Aspirin (Pain Relief)</option>
            <option value="caffeine">Caffeine (Stimulant)</option>
            <option value="custom">Custom (Manual Entry)</option>
          </select>
        </div>

        <div className="form-group">
          <label>Starting SMILES Seed (Manual)</label>
          <input
            type="text"
            value={config.seed_smiles}
            onChange={e => setConfig({ ...config, seed_smiles: e.target.value, preset: 'custom' })}
            style={{ fontFamily: 'monospace', fontSize: '0.75rem' }}
          />
        </div>

        <div className="form-group">
          <label>Generations</label>
          <input type="number" value={config.generations} onChange={e => setConfig({ ...config, generations: parseInt(e.target.value) })} />
        </div>

        <div className="form-group">
          <label>Island Pools</label>
          <input type="number" value={config.islands} onChange={e => setConfig({ ...config, islands: parseInt(e.target.value) })} />
        </div>

        <div className="form-group">
          <label>Population Size</label>
          <input type="number" value={config.pop_size} onChange={e => setConfig({ ...config, pop_size: parseInt(e.target.value) })} />
        </div>

        <button
          className="btn-primary"
          onClick={handleStart}
          disabled={evolving}
          style={{ marginTop: '1rem', opacity: evolving ? 0.5 : 1 }}
        >
          {evolving ? 'Evolving...' : 'Launch Generator'}
        </button>

        {currentBest && (
          <div style={{ marginTop: '2rem', padding: '1rem', background: 'rgba(0,0,0,0.3)', borderRadius: '8px' }}>
            <h3 style={{ fontSize: '0.9rem', color: 'var(--text-secondary)', marginBottom: '0.5rem' }}>Top Scaffold</h3>
            <p style={{ wordBreak: 'break-all', fontSize: '0.85rem' }}>{currentBest.smiles}</p>
          </div>
        )}
      </aside>

      <main className="main-content">
        <div className="stats-grid">
          <div className="stat-card glass-panel">
            <h3>Current Generation</h3>
            <div className="value">{generations.length}</div>
          </div>
          <div className="stat-card glass-panel">
            <h3>Best Fitness</h3>
            <div className="value" style={{ color: currentBest?.fitness > 0.4 ? 'var(--success-color)' : 'white' }}>
              {currentBest ? (currentBest.fitness * 100).toFixed(1) + '%' : '--'}
            </div>
          </div>
          <div className="stat-card glass-panel">
            <h3>Worst-Case Energy</h3>
            <div className="value">
              {currentBest
                ? Math.max(...Object.entries(currentBest.parts).filter(([k]) => k.startsWith('raw')).map(x => x[1])).toFixed(2)
                : '--'}
            </div>
          </div>
        </div>

        <div className="glass-panel" style={{ display: 'flex', flexDirection: 'column', padding: '1.5rem', height: '100%' }}>
          <h2 style={{ fontSize: '1.1rem', marginBottom: '1rem', display: 'flex', justifyContent: 'space-between' }}>
            Interactive 3D Visualization
            <span style={{ fontSize: '0.8rem', color: 'var(--text-secondary)' }}>Powered by 3Dmol.js</span>
          </h2>
          <div id="viewer3d" className="viewer-container" style={{ position: 'relative' }}>
            {!currentBest && !evolving && (
              <div style={{ color: 'var(--text-secondary)' }}>Start the evolution engine to view molecules.</div>
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
