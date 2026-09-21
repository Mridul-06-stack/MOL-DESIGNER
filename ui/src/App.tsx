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
  active_targets?: string[]
}

interface RedTeamAlert {
  type: 'redteam_alert'
  generation: number
  mutation: string
  mutation_short: string
  mutation_name?: string
  exon?: string
  mechanism?: string
  clinical_context?: string
  adaptation_guidance?: string
  trigger_smiles: string
  affinity_before: number
  active_targets: string[]
  message: string
}

// ── Sleek SVG Icons (No Emojis) ──────────────────────────────────────
const IconShield = () => (
  <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" style={{ flexShrink: 0 }}>
    <path d="M12 22s8-4 8-10V5l-8-3-8 3v7c0 6 8 10 8 10z"/>
  </svg>
)

const IconAlertTriangle = () => (
  <svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" style={{ flexShrink: 0 }}>
    <path d="m21.73 18-8-14a2 2 0 0 0-3.48 0l-8 14A2 2 0 0 0 4 21h16a2 2 0 0 0 1.73-3Z"/>
    <line x1="12" y1="9" x2="12" y2="13"/>
    <line x1="12" y1="17" x2="12.01" y2="17"/>
  </svg>
)

const IconDownload = () => (
  <svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" style={{ flexShrink: 0 }}>
    <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"/>
    <polyline points="7 10 12 15 17 10"/>
    <line x1="12" y1="15" x2="12" y2="3"/>
  </svg>
)

const IconFileSpreadsheet = () => (
  <svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" style={{ flexShrink: 0 }}>
    <path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"/>
    <polyline points="14 2 14 8 20 8"/>
    <line x1="8" y1="13" x2="16" y2="13"/>
    <line x1="8" y1="17" x2="16" y2="17"/>
    <line x1="10" y1="9" x2="8" y2="9"/>
  </svg>
)

const IconPlay = () => (
  <svg width="12" height="12" viewBox="0 0 24 24" fill="currentColor" style={{ flexShrink: 0 }}>
    <polygon points="5 3 19 12 5 21 5 3"/>
  </svg>
)

const CLINICAL_TARGET_META: Record<string, { label: string; badgeClass: string; indClass: string; role: string }> = {
  WT: { label: 'EGFR WT', badgeClass: 'target-pill-wt', indClass: 'indicator-wt', role: 'Native Kinase Domain' },
  L858R: { label: 'L858R', badgeClass: 'target-pill-driver', indClass: 'indicator-driver', role: 'Exon 21 Driver Mutation' },
  T790M: { label: 'T790M', badgeClass: 'target-pill-gatekeeper', indClass: 'indicator-gatekeeper', role: 'Gatekeeper Steric Clash' },
  C797S: { label: 'C797S', badgeClass: 'target-pill-covalent', indClass: 'indicator-covalent', role: 'Covalent Null Resistance' },
  L718Q: { label: 'L718Q', badgeClass: 'target-pill-ploop', indClass: 'indicator-ploop', role: 'P-Loop Polar Shift' },
  G724S: { label: 'G724S', badgeClass: 'target-pill-ploop', indClass: 'indicator-ploop', role: 'ATP Loop Conformation' },
}

function FitnessTrajectoryChart({ history }: { history: Array<{ generation: number; best: number; mean: number }> }) {
  if (history.length < 2) {
    return (
      <div style={{ height: '140px', display: 'flex', alignItems: 'center', justifyContent: 'center', color: 'var(--text-secondary)', fontSize: '0.85rem' }}>
        Run evolution to stream the fitness trajectory curve across generations.
      </div>
    )
  }

  const width = 600
  const height = 150
  const padding = { top: 15, right: 20, bottom: 25, left: 45 }
  const innerWidth = width - padding.left - padding.right
  const innerHeight = height - padding.top - padding.bottom

  const maxGen = Math.max(...history.map(h => h.generation), 1)
  const getX = (gen: number) => padding.left + (gen / maxGen) * innerWidth
  const getY = (val: number) => padding.top + innerHeight - Math.max(0, Math.min(val, 1)) * innerHeight

  const bestPoints = history.map(h => `${getX(h.generation)},${getY(h.best)}`).join(" ")
  const meanPoints = history.map(h => `${getX(h.generation)},${getY(h.mean)}`).join(" ")
  const areaPoints = `${getX(0)},${padding.top + innerHeight} ${bestPoints} ${getX(maxGen)},${padding.top + innerHeight}`

  const latest = history[history.length - 1]

  return (
    <div className="chart-container">
      <svg viewBox={`0 0 ${width} ${height}`} className="chart-svg" preserveAspectRatio="none">
        <defs>
          <linearGradient id="bestGradient" x1="0" y1="0" x2="0" y2="1">
            <stop offset="0%" stopColor="#38bdf8" stopOpacity="0.3" />
            <stop offset="100%" stopColor="#38bdf8" stopOpacity="0.0" />
          </linearGradient>
        </defs>

        {/* Grid lines */}
        {[0, 0.25, 0.5, 0.75, 1.0].map(val => (
          <g key={val}>
            <line
              x1={padding.left}
              y1={getY(val)}
              x2={width - padding.right}
              y2={getY(val)}
              stroke="rgba(255,255,255,0.07)"
              strokeDasharray="3,3"
            />
            <text
              x={padding.left - 8}
              y={getY(val) + 3}
              fill="#9ca3af"
              fontSize="9"
              textAnchor="end"
              fontFamily="sans-serif"
            >
              {(val * 100).toFixed(0)}%
            </text>
          </g>
        ))}

        {/* Area fill */}
        <polygon points={areaPoints} fill="url(#bestGradient)" />

        {/* Lines */}
        <polyline fill="none" stroke="#c084fc" strokeWidth="2" strokeDasharray="4,4" points={meanPoints} />
        <polyline fill="none" stroke="#38bdf8" strokeWidth="2.5" points={bestPoints} />

        {/* Current Dots */}
        <circle cx={getX(latest.generation)} cy={getY(latest.best)} r="4" fill="#38bdf8" stroke="#ffffff" strokeWidth="1.5" />
        <circle cx={getX(latest.generation)} cy={getY(latest.mean)} r="3.5" fill="#c084fc" stroke="#ffffff" strokeWidth="1.5" />

        {/* X Axis Labels */}
        <text x={padding.left} y={height - 6} fill="#9ca3af" fontSize="9" textAnchor="middle">Gen 0</text>
        <text x={width - padding.right} y={height - 6} fill="#9ca3af" fontSize="9" textAnchor="middle">Gen {maxGen}</text>
      </svg>
    </div>
  )
}

export default function App() {
  const [evolving, setEvolving] = useState(false)
  const [generations, setGenerations] = useState<GenerationEvent[]>([])
  const [currentBest, setCurrentBest] = useState<ScoredBest | null>(null)
  const [selectedCandidate, setSelectedCandidate] = useState<ScoredBest | null>(null)
  const [statusText, setStatusText] = useState<string>('')
  const [enableRedTeam, setEnableRedTeam] = useState(true)
  const [activeTargets, setActiveTargets] = useState<string[]>(["WT", "L858R"])
  const [latestAlert, setLatestAlert] = useState<RedTeamAlert | null>(null)
  const [topCandidates, setTopCandidates] = useState<ScoredBest[]>([])
  const [historyTimeline, setHistoryTimeline] = useState<Array<{ generation: number; best: number; mean: number }>>([])
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
    setSelectedCandidate(null)
    setLatestAlert(null)
    setTopCandidates([])
    setHistoryTimeline([])
    const initialTargets = enableRedTeam ? ["WT", "L858R"] : ["WT", "L858R", "T790M", "C797S"]
    setActiveTargets(initialTargets)
    setStatusText('Connecting to evolution engine...')

    try {
      // Stream generation events via SSE
      const res = await fetch('http://localhost:8000/api/evolve', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ 
          ...config, 
          targets: initialTargets,
          enable_redteam: enableRedTeam,
          redteam_threshold: -7.5,
          redteam_interval: 2,
          redteam_max_mutations: 3
        })
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
      setStatusText(enableRedTeam ? 'Simulating closed-loop adversarial evolution...' : 'Evolving chemistry...')

      while (true) {
        const { value, done } = await reader.read()
        if (done) break

        buffer += decoder.decode(value, { stream: true })
        const parts = buffer.split('\n\n')
        buffer = parts.pop() || ''

        for (const part of parts) {
          if (part.startsWith('data: ')) {
            try {
              const rawData = JSON.parse(part.slice(6))
              
              if (rawData.type === 'redteam_alert') {
                const alert: RedTeamAlert = rawData
                setLatestAlert(alert)
                if (alert.active_targets) {
                  setActiveTargets(alert.active_targets)
                }
              } else {
                const data: GenerationEvent = rawData
                setGenerations(prev => [...prev, data])
                if (data.active_targets) {
                  setActiveTargets(data.active_targets)
                }
                if (data.top) {
                  setTopCandidates(data.top.slice(0, 5))
                }
                if (data.history && data.best) {
                  setHistoryTimeline([
                    ...data.history,
                    { generation: data.generation, best: data.best.fitness, mean: data.mean_fitness }
                  ])
                }
                if (data.best) {
                  setCurrentBest(data.best)
                  setSelectedCandidate(data.best)
                  render3D(data.best.smiles)
                }
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

  const handleSelectCandidate = (candidate: ScoredBest) => {
    setSelectedCandidate(candidate)
    render3D(candidate.smiles)
  }

  const handleExportCSV = () => {
    if (topCandidates.length === 0 && !currentBest) return
    const listToExport = topCandidates.length > 0 ? topCandidates : [currentBest!]
    const headers = "Rank,SMILES,Fitness,QED,SA,MW,LogP,Worst_Energy_kcal_mol\n"
    const rows = listToExport.map((c, i) => {
      const rank = i + 1
      const smi = `"${c.smiles}"`
      const fit = (c.fitness * 100).toFixed(1) + "%"
      const qed = c.raw?.qed?.toFixed(4) || c.parts?.qed?.toFixed(4) || ""
      const sa = c.raw?.sa?.toFixed(2) || ""
      const mw = c.raw?.mw?.toFixed(2) || ""
      const logp = c.raw?.logp?.toFixed(2) || ""
      const energy = c.raw?.dock_worst != null ? c.raw.dock_worst.toFixed(2) : ""
      return `${rank},${smi},${fit},${qed},${sa},${mw},${logp},${energy}`
    }).join("\n")

    const blob = new Blob([headers + rows], { type: "text/csv;charset=utf-8;" })
    const url = URL.createObjectURL(blob)
    const link = document.createElement("a")
    link.href = url
    link.setAttribute("download", `alchemist_top_candidates_gen${generations.length}.csv`)
    document.body.appendChild(link)
    link.click()
    document.body.removeChild(link)
  }

  const handleDownloadSDF = async () => {
    const smi = selectedCandidate?.smiles || currentBest?.smiles
    if (!smi) return
    try {
      const res = await fetch(`http://localhost:8000/api/molblock?smiles=${encodeURIComponent(smi)}`)
      if (!res.ok) return
      const sdf = await res.text()
      const blob = new Blob([sdf], { type: "chemical/x-mdl-sdfile;charset=utf-8;" })
      const url = URL.createObjectURL(blob)
      const link = document.createElement("a")
      link.href = url
      link.setAttribute("download", `alchemist_candidate_conformer.sdf`)
      document.body.appendChild(link)
      link.click()
      document.body.removeChild(link)
    } catch (e) {
      console.error("Failed to download SDF:", e)
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

  const activeMolecule = selectedCandidate || currentBest
  const worstEnergy = activeMolecule?.raw?.dock_worst != null
    ? activeMolecule.raw.dock_worst.toFixed(2)
    : '--'

  return (
    <div className="dashboard-container">
      <header>
        <div className="header-brand">
          <div className="header-logo-row">
            <svg width="22" height="22" viewBox="0 0 24 24" fill="none" stroke="url(#headerGrad)" strokeWidth="2.2" strokeLinecap="round" strokeLinejoin="round">
              <defs>
                <linearGradient id="headerGrad" x1="0" y1="0" x2="1" y2="1">
                  <stop offset="0%" stopColor="#60a5fa" />
                  <stop offset="100%" stopColor="#c084fc" />
                </linearGradient>
              </defs>
              <circle cx="12" cy="12" r="3"/>
              <circle cx="19" cy="5" r="2"/>
              <circle cx="5" cy="19" r="2"/>
              <circle cx="5" cy="5" r="2"/>
              <circle cx="19" cy="19" r="2"/>
              <line x1="12" y1="9" x2="19" y2="5"/>
              <line x1="12" y1="15" x2="5" y2="19"/>
              <line x1="9.5" y1="10" x2="5" y2="5"/>
              <line x1="14.5" y1="14" x2="19" y2="19"/>
            </svg>
            <h1>Alchemist</h1>
            <span className="version-tag">v0.1.0</span>
          </div>
          <p className="header-subtitle">
            Closed-Loop Adversarial AI Drug Discovery &bull; EGFR Resistance Overcoming
          </p>
        </div>
        {evolving && (
          <div className={`header-status-pill ${enableRedTeam ? 'status-adversarial' : 'status-standard'}`}>
            <span className="live-beacon" />
            <span>{enableRedTeam ? 'Adversarial Co-Evolution Loop' : 'Generative Evolution'}</span>
          </div>
        )}
      </header>

      <aside className="sidebar glass-panel">
        <h2 className="sidebar-section-title">
          Simulation Parameters
        </h2>

        {/* Adversarial Red-Team Toggle */}
        <div className="toggle-wrapper">
          <div className="toggle-label">
            <span className="toggle-title">
              <IconShield /> Red-Team Mode
            </span>
            <span className="toggle-subtitle">Simulate tumor resistance</span>
          </div>
          <label className="toggle-switch">
            <input
              type="checkbox"
              checked={enableRedTeam}
              onChange={e => setEnableRedTeam(e.target.checked)}
              disabled={evolving}
            />
            <span className="toggle-slider"></span>
          </label>
        </div>

        {/* Active Targets Panel */}
        <div className="form-group" style={{ marginBottom: '1.25rem' }}>
          <label style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
            <span>Active Target Panel</span>
            <span style={{ fontSize: '0.75rem', color: 'var(--text-secondary)' }}>{activeTargets.length} variants</span>
          </label>
          <div className="targets-container">
            {activeTargets.map(tgt => {
              const meta = CLINICAL_TARGET_META[tgt]
              const isResist = tgt.startsWith('RESIST_')
              let label = meta ? meta.label : tgt
              let badgeClass = meta ? meta.badgeClass : (isResist ? 'target-pill-resist' : tgt !== 'WT' ? 'target-pill-driver' : 'target-pill-wt')
              let indClass = meta ? meta.indClass : (isResist ? 'indicator-resist' : tgt !== 'WT' ? 'indicator-driver' : 'indicator-wt')
              let tooltip = meta ? `${tgt}: ${meta.role}` : tgt

              if (isResist) {
                const parts = tgt.split('_')
                label = parts[2] ? `RESIST-${parts[2]}` : tgt.slice(0, 10)
                tooltip = `Synthetic escape mutation: ${tgt}`
              }
              return (
                <span
                  key={tgt}
                  className={`target-pill ${badgeClass}`}
                  title={tooltip}
                >
                  <span className={`target-indicator ${indClass}`} />
                  {label}
                </span>
              )
            })}
          </div>
        </div>

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
          id="launch-evolution-btn"
          className="btn-primary"
          onClick={handleStart}
          disabled={evolving}
          style={{ 
            marginTop: '1rem', 
            opacity: evolving ? 0.7 : 1,
            background: enableRedTeam ? 'linear-gradient(135deg, #e11d48, #7c3aed)' : undefined
          }}
        >
          {evolving ? (
            <span className="btn-inner">
              <span className="btn-spinner" />
              <span>{enableRedTeam ? 'Simulating Adversarial Loop...' : 'Evolving Chemistry...'}</span>
            </span>
          ) : (
            <span className="btn-inner">
              {enableRedTeam ? <IconShield /> : <IconPlay />}
              <span>{enableRedTeam ? 'Launch Adversarial Engine' : 'Launch Generator'}</span>
            </span>
          )}
        </button>

        {statusText && (
          <p style={{ fontSize: '0.8rem', color: 'var(--text-secondary)', marginTop: '0.5rem', textAlign: 'center' }}>
            {statusText}
          </p>
        )}

        {activeMolecule && (
          <div style={{ marginTop: '1.5rem', padding: '1rem', background: 'rgba(0,0,0,0.3)', borderRadius: '8px', border: '1px solid rgba(255,255,255,0.08)' }}>
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '0.5rem' }}>
              <h3 style={{ fontSize: '0.85rem', color: 'var(--text-secondary)', textTransform: 'uppercase', letterSpacing: '0.05em' }}>
                {selectedCandidate && selectedCandidate !== currentBest ? 'Inspecting Candidate' : 'Current Best Lead'}
              </h3>
              <span className="pill" style={{ margin: 0, fontSize: '0.7rem', background: 'rgba(59,130,246,0.2)', color: 'var(--accent-blue)' }}>
                {(activeMolecule.fitness * 100).toFixed(1)}% Fit
              </span>
            </div>
            <p style={{ wordBreak: 'break-all', fontSize: '0.8rem', fontFamily: 'monospace', color: '#60a5fa' }}>
              {activeMolecule.smiles}
            </p>
            {activeMolecule.raw && (
              <div style={{ marginTop: '0.75rem', fontSize: '0.75rem', color: 'var(--text-secondary)', display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '0.25rem' }}>
                <div>MW: {activeMolecule.raw.mw?.toFixed(1)} Da</div>
                <div>LogP: {activeMolecule.raw.logp?.toFixed(2)}</div>
                <div>QED: {activeMolecule.raw.qed?.toFixed(2)}</div>
                <div>SA: {activeMolecule.raw.sa?.toFixed(2)}</div>
              </div>
            )}

            {/* Per-Variant Affinity Breakdown */}
            {activeMolecule.raw?.dock && Object.keys(activeMolecule.raw.dock).length > 0 && (
              <div className="variant-breakdown">
                <div className="variant-breakdown-title">Variant Affinity Profile (kcal/mol)</div>
                <div className="variant-energy-list">
                  {Object.entries(activeMolecule.raw.dock).map(([variant, energy]) => {
                    const meta = CLINICAL_TARGET_META[variant]
                    const isResist = variant.startsWith('RESIST_')
                    const displayName = meta ? meta.label : (isResist ? (variant.split('_')[2] ? `RES-${variant.split('_')[2]}` : variant.slice(0, 10)) : variant)
                    const role = meta ? meta.role : (isResist ? 'Synthetic Escape' : 'Kinase Variant')
                    const clamped = Math.min(Math.max(-energy, 0), 10)
                    const pct = (clamped / 10) * 100
                    const barColor = energy < -8 ? 'var(--success-color)' : energy < -6 ? '#f59e0b' : '#ef4444'
                    
                    let textColor = '#d1d5db'
                    if (variant === 'T790M') textColor = '#fbbf24'
                    else if (variant === 'C797S') textColor = '#f87171'
                    else if (variant === 'L718Q' || variant === 'G724S') textColor = '#c084fc'
                    else if (variant === 'L858R') textColor = '#38bdf8'
                    else if (variant === 'WT') textColor = '#34d399'
                    else if (isResist) textColor = '#fca5a5'

                    return (
                      <div key={variant} className="variant-energy-row">
                        <span 
                          style={{ 
                            width: '82px', 
                            overflow: 'hidden', 
                            textOverflow: 'ellipsis', 
                            whiteSpace: 'nowrap', 
                            color: textColor,
                            fontWeight: 500
                          }} 
                          title={`${variant}: ${role}`}
                        >
                          {displayName}
                        </span>
                        <div className="variant-energy-bar-wrap">
                          <div className="variant-energy-bar" style={{ width: `${pct}%`, backgroundColor: barColor }} />
                        </div>
                        <span style={{ fontFamily: 'monospace', width: '55px', textAlign: 'right', color: barColor }}>
                          {energy.toFixed(2)}
                        </span>
                      </div>
                    )
                  })}
                </div>
              </div>
            )}
          </div>
        )}
      </aside>

      <main className="main-content">
        {/* Red-Team Alert Banner */}
        {latestAlert && (
          <div className="redteam-banner">
            <div className="redteam-header">
              <div className="redteam-badge">
                <IconAlertTriangle />
                <span>Clinical Mutational Escape</span>
              </div>
              <span className="redteam-gen-tag">
                Generation {latestAlert.generation} • {latestAlert.exon || 'Exon 20'}
              </span>
            </div>
            <div className="redteam-title">
              {latestAlert.mutation_name || `Target Escape: ${latestAlert.mutation_short}`}
            </div>
            <div className="redteam-desc">
              {latestAlert.clinical_context && (
                <div style={{ marginBottom: '0.4rem' }}>
                  <strong style={{ color: '#f1f5f9' }}>Oncology Context:</strong> {latestAlert.clinical_context}
                </div>
              )}
              <div style={{ marginBottom: '0.4rem' }}>
                <strong style={{ color: '#f1f5f9' }}>Biophysical Mechanism:</strong> {latestAlert.mechanism || latestAlert.message}
              </div>
              {latestAlert.adaptation_guidance && (
                <div className="redteam-guidance">
                  <span className="redteam-guidance-badge">Generator Action:</span> {latestAlert.adaptation_guidance}
                </div>
              )}
            </div>
            <div className="redteam-footer">
              <div className="redteam-mutation-tag">
                Ensemble Injected: {latestAlert.mutation}
              </div>
              <div className="redteam-trigger-tag">
                Trigger Lead Affinity: {latestAlert.affinity_before.toFixed(2)} kcal/mol
              </div>
            </div>
          </div>
        )}

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
            <h3>Worst-Case Energy</h3>
            <div className="value" style={{ color: activeMolecule?.raw?.dock_worst && activeMolecule.raw.dock_worst < -8 ? 'var(--success-color)' : 'white' }}>
              {worstEnergy} <span style={{ fontSize: '1rem', color: 'var(--text-secondary)' }}>kcal/mol</span>
            </div>
          </div>
        </div>

        {/* Generational Fitness Trajectory Chart */}
        <div className="glass-panel" style={{ padding: '1.25rem' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
            <h2 style={{ fontSize: '1rem', fontWeight: 600 }}>Generational Fitness Trajectory</h2>
            <div className="chart-legend">
              <div className="legend-item">
                <span className="legend-dot" style={{ background: '#38bdf8' }} />
                <span>Best Fitness</span>
              </div>
              <div className="legend-item">
                <span className="legend-dot" style={{ background: '#c084fc' }} />
                <span>Mean Population</span>
              </div>
            </div>
          </div>
          <FitnessTrajectoryChart history={historyTimeline} />
        </div>

        {/* 3D Molecular Conformer Viewer */}
        <div className="glass-panel" style={{ display: 'flex', flexDirection: 'column', padding: '1.5rem', minHeight: '500px' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
            <div>
              <h2 style={{ fontSize: '1.1rem', fontWeight: 600 }}>
                Interactive 3D Molecular Conformer
              </h2>
              <span style={{ fontSize: '0.8rem', color: 'var(--text-secondary)' }}>
                {selectedCandidate && selectedCandidate !== currentBest ? 'Displaying Selected Candidate from Leaderboard' : 'Displaying Lead Candidate'}
              </span>
            </div>
            {activeMolecule && (
              <div className="export-toolbar" style={{ margin: 0 }}>
                <button className="btn-secondary" onClick={handleDownloadSDF}>
                  <IconDownload />
                  <span>Download 3D (.sdf)</span>
                </button>
                <button className="btn-secondary" onClick={handleExportCSV}>
                  <IconFileSpreadsheet />
                  <span>Export Candidates (.csv)</span>
                </button>
              </div>
            )}
          </div>

          <div id="viewer3d" className="viewer-container" style={{ position: 'relative', width: '100%', height: '420px' }}>
            {!currentBest && !evolving && (
              <div style={{ color: 'var(--text-secondary)', textAlign: 'center', padding: '2rem' }}>
                Click <strong>"Launch Generator"</strong> or <strong>"Launch Adversarial Engine"</strong> to begin evolution and render 3D conformers.
              </div>
            )}
            {evolving && !currentBest && (
              <div className="loading-spinner"></div>
            )}
          </div>
        </div>

        {/* Interactive Top Candidates Leaderboard */}
        {topCandidates.length > 0 && (
          <div className="glass-panel" style={{ padding: '1.25rem' }}>
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
              <div>
                <h2 style={{ fontSize: '1rem', fontWeight: 600 }}>Top Surviving Candidates (Click to View 3D)</h2>
                <p style={{ fontSize: '0.75rem', color: 'var(--text-secondary)' }}>
                  Select any molecular candidate to switch the 3D conformer viewer above.
                </p>
              </div>
              <button className="btn-secondary" onClick={handleExportCSV}>
                <IconFileSpreadsheet />
                <span>Export All (CSV)</span>
              </button>
            </div>

            <div className="leaderboard-container">
              {topCandidates.map((candidate, idx) => {
                const isSelected = selectedCandidate?.smiles === candidate.smiles
                const worst = candidate.raw?.dock_worst != null ? candidate.raw.dock_worst.toFixed(2) : '--'
                return (
                  <div
                    key={candidate.smiles}
                    className={`leaderboard-row ${isSelected ? 'active' : ''}`}
                    onClick={() => handleSelectCandidate(candidate)}
                  >
                    <div className="leaderboard-left">
                      <div className={`rank-badge ${idx === 0 ? 'rank-1' : idx === 1 ? 'rank-2' : idx === 2 ? 'rank-3' : ''}`}>
                        #{idx + 1}
                      </div>
                      <div className="leaderboard-smiles" title={candidate.smiles}>
                        {candidate.smiles}
                      </div>
                    </div>

                    <div className="leaderboard-right">
                      <div className="metric-pill">
                        <span className="metric-pill-label">Fitness</span>
                        <span className="metric-pill-val" style={{ color: candidate.fitness > 0.6 ? 'var(--success-color)' : 'white' }}>
                          {(candidate.fitness * 100).toFixed(1)}%
                        </span>
                      </div>
                      <div className="metric-pill">
                        <span className="metric-pill-label">Worst Energy</span>
                        <span className="metric-pill-val" style={{ color: candidate.raw?.dock_worst && candidate.raw.dock_worst < -8 ? 'var(--success-color)' : '#f59e0b' }}>
                          {worst}
                        </span>
                      </div>
                      <div className="metric-pill">
                        <span className="metric-pill-label">QED / SA</span>
                        <span className="metric-pill-val" style={{ color: 'var(--text-secondary)' }}>
                          {candidate.raw?.qed?.toFixed(2) || '--'} / {candidate.raw?.sa?.toFixed(1) || '--'}
                        </span>
                      </div>
                    </div>
                  </div>
                )
              })}
            </div>
          </div>
        )}
      </main>
    </div>
  )
}



