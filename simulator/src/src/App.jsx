import { useState } from 'react'
import FertilityModelFit   from './components/FertilityModelFit.jsx'
import GenerationParams    from './components/GenerationParams.jsx'
import KinTypeSelector     from './components/KinTypeSelector.jsx'
import SimulationControls  from './components/SimulationControls.jsx'
import ResultsChart        from './components/ResultsChart.jsx'
import SummaryStats        from './components/SummaryStats.jsx'
import HistoricalTrends    from './components/HistoricalTrends.jsx'
import CohortValidationTable from './components/CohortValidationTable.jsx'
import { simulateKin }     from './lib/kinSimulator.js'

// Default params per model type
const DEFAULTS = {
  zinb: {
    focal:       { mu: 3.213, theta: 19.536, pi0: 0.056 },  // 1990 cohort
    parent:      { mu: 2.530, theta:  3.652, pi0: 0.066 },  // 1970 cohort
    grandparent: { mu: 2.943, theta:  2.372, pi0: 0.043 },  // 1950 cohort
  },
  poisson: {
    focal:       { mu: 3.034 },
    parent:      { mu: 2.363 },
    grandparent: { mu: 2.818 },
  },
  fixed: {
    focal:       { fertMean: 3.034 },
    parent:      { fertMean: 2.363, sibMean: 3.314 },
    grandparent: { sibMean:  4.192 },
  },
}

const DEFAULT_KIN = {
  children: true, siblings: true, auntsUncles: true,
  cousins: true, niecesNephews: false,
}

const TABS = [
  { id: 'fit', label: 'Fertility Fit',
    desc: 'How well does each model capture the empirical fertility and sibling distributions?' },
  { id: 'kin', label: 'Kin Counts',
    desc: 'Simulate kin count distributions under a chosen fertility model across three generations.' },
]

const MODEL_OPTIONS = [
  { key: 'zinb',    label: 'ZINB',    desc: 'Zero-inflated negative binomial — captures overdispersion and excess zeros' },
  { key: 'poisson', label: 'Poisson', desc: 'Equidispersed — variance equals mean' },
  { key: 'fixed',   label: 'Fixed',   desc: 'Fixed means — fertility and sibling sizes set to their empirical means' },
]

export default function App() {
  const [activeTab,   setActiveTab]   = useState('fit')
  const [model,       setModel]       = useState('zinb')
  const [genParams,   setGenParams]   = useState(DEFAULTS.zinb)
  const [selectedKin, setSelectedKin] = useState(DEFAULT_KIN)
  const [seed,        setSeed]        = useState(42)
  const [results,     setResults]     = useState(null)
  const [isRunning,   setIsRunning]   = useState(false)
  const [runInfo,     setRunInfo]     = useState(null)
  const [nSim,        setNSim]        = useState(100000)
  const [selectedYear, setSelectedYear] = useState(1990)

  function handleModelChange(m) {
    setModel(m)
    setGenParams(DEFAULTS[m])
    setResults(null)
  }

  function handleGenChange(genKey, params) {
    setGenParams(prev => ({ ...prev, [genKey]: params }))
    setResults(null)
  }

  function handleRun() {
    setIsRunning(true)
    setTimeout(() => {
      const t0  = performance.now()
      const res = simulateKin({
        model,
        focal:       genParams.focal,
        parent:      genParams.parent,
        grandparent: genParams.grandparent,
        nSim,
        seed,
      })
      const ms = (performance.now() - t0).toFixed(0)
      setResults(res)
      setRunInfo({ ms, nSim, model })
      setIsRunning(false)
    }, 10)
  }

  const activeTabMeta = TABS.find(t => t.id === activeTab)

  return (
    <div className="app">
      <header className="app-header">
        <div className="app-header-text">
          <h1>Kin Counts Simulator</h1>
          <p>Bridging fertility and sibling distributions — interactive companion to the paper</p>
        </div>
      </header>

      <nav className="tab-nav">
        {TABS.map(tab => (
          <button
            key={tab.id}
            className={`tab-btn ${activeTab === tab.id ? 'active' : ''}`}
            onClick={() => setActiveTab(tab.id)}
          >
            {tab.label}
          </button>
        ))}
      </nav>
      <p className="tab-desc">{activeTabMeta.desc}</p>

      {/* ══ Fertility Fit ══ */}
      {activeTab === 'fit' && (
        <div className="tab-content">
          <FertilityModelFit selectedYear={selectedYear} onYearChange={setSelectedYear} />
        </div>
      )}

      {/* ══ Kin Counts ══ */}
      {activeTab === 'kin' && (
        <div className="tab-content">

          <HistoricalTrends />
          <div className="kin-divider" />
          <CohortValidationTable />
          <div className="kin-divider" />

          <div className="app-body">
            <aside className="controls-panel">

              {/* Step 1: Model selection */}
              <div className="control-group">
                <span className="control-label">Fertility Model</span>
                {MODEL_OPTIONS.map(({ key, label, desc }) => (
                  <label key={key} className="radio-label" style={{ alignItems: 'flex-start', gap: '0.5rem' }}>
                    <input
                      type="radio"
                      name="model"
                      value={key}
                      checked={model === key}
                      onChange={() => handleModelChange(key)}
                      style={{ marginTop: '0.2rem' }}
                    />
                    <span>
                      <strong>{label}</strong>
                      <span className="model-hint" style={{ display: 'block' }}>{desc}</span>
                    </span>
                  </label>
                ))}
              </div>

              {/* Step 2: Per-generation parameters */}
              <div className="control-group">
                <span className="control-label">Generation Parameters</span>
                {['focal', 'parent', 'grandparent'].map(genKey => (
                  <GenerationParams
                    key={`${model}-${genKey}`}
                    genKey={genKey}
                    model={model}
                    params={genParams[genKey]}
                    onChange={p => handleGenChange(genKey, p)}
                  />
                ))}
              </div>

              {/* Step 3: Kin types + run */}
              <KinTypeSelector selectedKin={selectedKin} onChange={setSelectedKin} />

              <SimulationControls
                onRun={handleRun}
                isRunning={isRunning}
                seed={seed}
                onSeedChange={setSeed}
                nSim={nSim}
                onNSimChange={v => { setNSim(v); setResults(null) }}
              />
            </aside>

            <main className="results-panel">
              {runInfo && !isRunning && (
                <p className="run-meta">
                  {runInfo.nSim.toLocaleString()} draws · {runInfo.model.toUpperCase()} · {runInfo.ms} ms
                </p>
              )}
              <h2>Kin Count Distributions</h2>
              <ResultsChart results={results} selectedKin={selectedKin} />
              <h2>Summary Statistics</h2>
              <SummaryStats results={results} selectedKin={selectedKin} />
              {!results && (
                <p className="results-hint">
                  Select a model, set parameters for each generation, then click <strong>Run Simulation</strong>.
                </p>
              )}
            </main>
          </div>
        </div>
      )}
    </div>
  )
}
