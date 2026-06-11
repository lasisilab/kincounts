const MODELS = [
  {
    value: 'poisson',
    label: 'Poisson',
    hint: 'Pure Poisson fertility. Variance = mean — severely underestimates overdispersion.',
  },
  {
    value: 'zinb',
    label: 'ZINB',
    hint: 'Zero-inflated negative binomial. θ controls overdispersion, π₀ excess zeros. Best fit to census data.',
  },
  {
    value: 'fixed',
    label: 'Fixed',
    hint: 'Deterministic mean only — every family has exactly μ children. Zero distributional variance.',
  },
]

export default function FertilityModelSelector({ model, onChange }) {
  const hint = MODELS.find(m => m.value === model)?.hint

  return (
    <div className="control-group">
      <span className="control-label">Fertility Model</span>
      <div className="radio-group">
        {MODELS.map(({ value, label }) => (
          <label key={value} className="radio-label">
            <input
              type="radio"
              name="model"
              value={value}
              checked={model === value}
              onChange={() => onChange(value)}
            />
            {label}
          </label>
        ))}
      </div>
      {hint && <p className="model-hint">{hint}</p>}
    </div>
  )
}
