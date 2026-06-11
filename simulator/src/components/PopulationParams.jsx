import { IPUMS_COHORTS } from '../lib/empiricalData.js'

function Slider({ label, value, min, max, step, onChange, display }) {
  const formatted = display ? display(value) : value
  return (
    <div className="slider-row">
      <div className="slider-header">
        <span className="slider-label">{label}</span>
        <span className="slider-value">{formatted}</span>
      </div>
      <input
        type="range"
        min={min}
        max={max}
        step={step}
        value={value}
        onChange={e => onChange(Number(e.target.value))}
      />
    </div>
  )
}

export default function PopulationParams({ params, onChange, onCohortSelect }) {
  const set = (key, val) => onChange({ ...params, [key]: val })

  function handlePreset(e) {
    const year = Number(e.target.value)
    if (!year) return
    const cohort = IPUMS_COHORTS.find(c => c.year === year)
    if (!cohort) return
    onChange({ ...params, mu: cohort.mu, theta: cohort.theta, pi0: cohort.pi0 })
    onCohortSelect(cohort)
  }

  const effectiveMean = ((1 - params.pi0) * params.mu).toFixed(3)

  return (
    <div className="control-group">
      <span className="control-label">ZINB Parameters</span>

      <div className="preset-row">
        <label className="preset-label">IPUMS Cohort Preset</label>
        <select className="preset-select" defaultValue="" onChange={handlePreset}>
          <option value="">— manual —</option>
          {IPUMS_COHORTS.map(c => (
            <option key={c.year} value={c.year}>
              {c.year} ({c.cohort})
            </option>
          ))}
        </select>
      </div>

      <Slider
        label="μ — NB component mean"
        value={params.mu}
        min={1} max={6} step={0.05}
        onChange={v => set('mu', v)}
        display={v => v.toFixed(2)}
      />
      <Slider
        label="θ — dispersion (higher = less dispersed)"
        value={params.theta}
        min={1} max={25} step={0.1}
        onChange={v => set('theta', v)}
        display={v => v.toFixed(1)}
      />
      <Slider
        label="π₀ — zero-inflation prob"
        value={params.pi0}
        min={0} max={0.2} step={0.01}
        onChange={v => set('pi0', v)}
        display={v => v.toFixed(2)}
      />
      <p className="effective-mean">
        Effective mean = (1−π₀)×μ = <strong>{effectiveMean}</strong>
        <br />Poisson and Fixed use the same mean for a fair comparison.
      </p>

      <Slider
        label="Simulations"
        value={params.nSim}
        min={1000} max={50000} step={1000}
        onChange={v => set('nSim', v)}
        display={v => v.toLocaleString()}
      />
    </div>
  )
}
