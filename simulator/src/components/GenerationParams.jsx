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

function NumberInput({ label, value, min, max, step, onChange }) {
  return (
    <div className="slider-row">
      <div className="slider-header">
        <span className="slider-label">{label}</span>
        <input
          className="gen-number-input"
          type="number" min={min} max={max} step={step}
          value={value}
          onChange={e => {
            const v = Number(e.target.value)
            if (!isNaN(v)) onChange(v)
          }}
        />
      </div>
    </div>
  )
}

const GEN_LABELS = {
  focal:       'Focal Generation',
  parent:      'Parent Generation',
  grandparent: 'Grandparent Generation',
}

const GEN_NOTES = {
  focal:       'X₀ — used for children and nieces/nephews',
  parent:      'X₁ — used for siblings and cousins',
  grandparent: 'X₂ — used for aunts and uncles',
}

export default function GenerationParams({ genKey, model, params, seed, onChange }) {
  const set = (key, val) => onChange({ ...params, [key]: val })

  return (
    <div className="gen-block">
      <div className="gen-block-header">
        <span className="gen-block-title">{GEN_LABELS[genKey]}</span>
        <span className="gen-block-note">{GEN_NOTES[genKey]}</span>
        {seed && (
          <span className="gen-block-seed">
            from {seed.year}{seed.cohort ? ` · b. ${seed.cohort}` : ''}
          </span>
        )}
      </div>

      {model === 'fixed' ? (
        <>
          {/* Fixed: independent fertility mean and sibling mean */}
          {(genKey === 'focal' || genKey === 'parent') && (
            <NumberInput
              label="Fertility mean"
              value={params.fertMean ?? 3.0}
              min={0.5} max={8} step={0.01}
              onChange={v => set('fertMean', v)}
            />
          )}
          {(genKey === 'parent' || genKey === 'grandparent') && (
            <NumberInput
              label="Sibling mean"
              value={params.sibMean ?? 3.5}
              min={0.5} max={10} step={0.01}
              onChange={v => set('sibMean', v)}
            />
          )}
        </>
      ) : model === 'poisson' ? (
        <Slider
          label="μ — mean fertility"
          value={params.mu ?? 3.0}
          min={0.5} max={6} step={0.05}
          onChange={v => set('mu', v)}
          display={v => v.toFixed(2)}
        />
      ) : (
        /* ZINB */
        <>
          <Slider
            label="μ — NB component mean"
            value={params.mu ?? 3.0}
            min={0.5} max={6} step={0.05}
            onChange={v => set('mu', v)}
            display={v => v.toFixed(2)}
          />
          <Slider
            label="θ — dispersion"
            value={params.theta ?? 5}
            min={0.5} max={25} step={0.1}
            onChange={v => set('theta', v)}
            display={v => v.toFixed(1)}
          />
          <Slider
            label="π₀ — zero-inflation"
            value={params.pi0 ?? 0.05}
            min={0} max={0.3} step={0.01}
            onChange={v => set('pi0', v)}
            display={v => v.toFixed(2)}
          />
          <p className="effective-mean">
            Effective mean = (1−π₀)×μ = {((1 - (params.pi0 ?? 0.05)) * (params.mu ?? 3.0)).toFixed(3)}
          </p>
        </>
      )}
    </div>
  )
}
