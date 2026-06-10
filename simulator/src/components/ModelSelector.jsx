const MODELS = [
  { key: 'zinb',    label: 'ZINB',    color: '#e8704a', note: 'overdispersed' },
  { key: 'poisson', label: 'Poisson', color: '#4f86c6', note: 'equidispersed' },
  { key: 'fixed',   label: 'Fixed',   color: '#8b8b8b', note: 'fixed means'  },
]

export default function ModelSelector({ selectedModels, onChange }) {
  function toggle(key) {
    // Prevent deselecting all
    const next = { ...selectedModels, [key]: !selectedModels[key] }
    if (!Object.values(next).some(Boolean)) return
    onChange(next)
  }

  return (
    <div className="control-group">
      <span className="control-label">Models to Compare</span>
      {MODELS.map(({ key, label, color, note }) => (
        <label key={key} className="checkbox-label">
          <input
            type="checkbox"
            checked={!!selectedModels[key]}
            onChange={() => toggle(key)}
          />
          <span style={{ color, fontWeight: 600 }}>{label}</span>
          <span className="model-check-note">{note}</span>
        </label>
      ))}
    </div>
  )
}
