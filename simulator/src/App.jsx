import { useState, useMemo } from 'react'
import FertilityModelFit     from './components/FertilityModelFit.jsx'
import GenerationParams      from './components/GenerationParams.jsx'
import KinTypeSelector       from './components/KinTypeSelector.jsx'
import SimulationControls    from './components/SimulationControls.jsx'
import ResultsChart          from './components/ResultsChart.jsx'
import SummaryStats          from './components/SummaryStats.jsx'
import HistoricalTrends      from './components/HistoricalTrends.jsx'
import CohortValidationTable from './components/CohortValidationTable.jsx'
import DataImport         from './components/DataImport.jsx'
import { simulateKin }       from './lib/kinSimulator.js'
import { IPUMS_DATASET, defaultsFromDataset } from './lib/datasets.js'

const DEFAULT_KIN = {
  children: true, siblings: true, auntsUncles: true,
  cousins: true, niecesNephews: false,
}

const TABS = [
  { id: 'fit', label: 'Fertility Fit',
    desc: 'The empirical fertility and sibling distributions from the data, and how the ZINB, Poisson, and fixed-mean models capture them.' },
  { id: 'sim', label: 'Simulator',
    desc: 'Simulate kin-count distributions under a chosen fertility model across three generations.' },
]

const MODEL_OPTIONS = [
  { key: 'zinb',    label: 'ZINB',    desc: 'Zero-inflated negative binomial — captures overdispersion and excess zeros' },
  { key: 'poisson', label: 'Poisson', desc: 'Equidispersed — variance equals mean' },
  { key: 'fixed',   label: 'Fixed',   desc: 'Fixed means — fertility and sibling sizes set to their empirical means' },
]

function newestYear(dataset) {
  return [...dataset.cohorts].sort((a, b) => a.year - b.year).at(-1).year
}

export default function App() {
  // The ACTIVE dataset. Defaults to IPUMS; an imported dataset replaces it
  // everywhere (Fertility Fit + Simulator) with no special-casing.
  const [dataset,     setDataset]     = useState(IPUMS_DATASET)
  const defaults = useMemo(() => defaultsFromDataset(dataset), [dataset])

  const [importOpen,  setImportOpen]  = useState(false)
  const [activeTab,   setActiveTab]   = useState('fit')
  const [model,       setModel]       = useState('zinb')
  const [genParams,   setGenParams]   = useState(() => defaultsFromDataset(IPUMS_DATASET).byModel.zinb)
  const [selectedKin, setSelectedKin] = useState(DEFAULT_KIN)
  const [seed,        setSeed]        = useState(42)
  const [results,     setResults]     = useState(null)
  const [isRunning,   setIsRunning]   = useState(false)
  const [runInfo,     setRunInfo]     = useState(null)
  const [nSim,        setNSim]        = useState(100000)
  const [selectedYear, setSelectedYear] = useState(() => newestYear(IPUMS_DATASET))

  // Swap the active dataset (IPUMS or imported) and reset everything that
  // depends on it: model defaults, selected year, and any prior simulation.
  function applyDataset(ds) {
    const d = defaultsFromDataset(ds)
    setDataset(ds)
    setModel('zinb')
    setGenParams(d.byModel.zinb)
    setSelectedYear(newestYear(ds))
    setResults(null)
    setImportOpen(false)
  }

  function handleModelChange(m) {
    setModel(m)
    setGenParams(defaults.byModel[m])
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
  const activeCohort  = dataset.cohorts.find(c => c.year === selectedYear)

  return (
    <div className="app">
      <header className="app-header">
        <div className="app-header-text">
          <h1>Kincounts</h1>
          <p>Fertility and kin-count distributions — interactive companion to the paper</p>
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

      <div className="dataset-bar">
        <span className="dataset-bar-label">Data source</span>
        <span className="dataset-bar-name">{dataset.label}</span>
        <button className="dataset-bar-btn" onClick={() => setImportOpen(o => !o)}>
          {importOpen ? 'Close importer' : 'Use your own data'}
        </button>
        {dataset.id !== 'ipums' && (
          <button className="dataset-bar-reset" onClick={() => applyDataset(IPUMS_DATASET)}>
            Reset to IPUMS
          </button>
        )}
      </div>

      {importOpen && (
        <DataImport
          baseDataset={IPUMS_DATASET}
          onLoad={applyDataset}
          onClose={() => setImportOpen(false)}
        />
      )}

      {/* ══ Fertility Fit ══ */}
      {activeTab === 'fit' && (
        <div className="tab-content">
          <FertilityModelFit
            dataset={dataset}
            selectedYear={selectedYear}
            onYearChange={setSelectedYear}
          />
          <div className="kin-divider" />
          <HistoricalTrends dataset={dataset} activeCohort={activeCohort} />
          <div className="kin-divider" />
          <CohortValidationTable dataset={dataset} />
        </div>
      )}

      {/* ══ Simulator ══ */}
      {activeTab === 'sim' && (
        <div className="tab-content">
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
                <p className="control-sublabel">
                  Defaults seeded from {dataset.label}: focal = {defaults.seeds.focal.year}, parent = {defaults.seeds.parent.year}, grandparent = {defaults.seeds.grandparent.year}.
                </p>
                {['focal', 'parent', 'grandparent'].map(genKey => (
                  <GenerationParams
                    key={`${dataset.id}-${model}-${genKey}`}
                    genKey={genKey}
                    model={model}
                    params={genParams[genKey]}
                    seed={defaults.seeds[genKey]}
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
