function Slider({ label, value, min, max, step, onChange, display }) {
  return (
    <div className="slider-row">
      <div className="slider-header">
        <span className="slider-label">{label}</span>
        <span className="slider-value">{display ? display(value) : value}</span>
      </div>
      <input
        type="range" min={min} max={max} step={step} value={value}
        onChange={e => onChange(Number(e.target.value))}
      />
    </div>
  )
}

export default function SimulationControls({ onRun, isRunning, seed, onSeedChange, nSim, onNSimChange }) {
  return (
    <div className="control-group">
      <span className="control-label">Run</span>
      <Slider
        label="Simulations"
        value={nSim}
        min={10000} max={200000} step={10000}
        onChange={onNSimChange}
        display={v => v.toLocaleString()}
      />
      <div className="seed-row">
        <label htmlFor="seed-input">Seed</label>
        <input
          id="seed-input"
          type="number"
          value={seed}
          min={0}
          onChange={e => onSeedChange(Number(e.target.value))}
        />
      </div>
      <button className="run-button" onClick={onRun} disabled={isRunning}>
        {isRunning ? 'Running…' : 'Run Simulation'}
      </button>
    </div>
  )
}
