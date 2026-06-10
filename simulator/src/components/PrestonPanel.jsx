// Displays Preston's identity: E[Y] = X̄ + σ²_X/X̄ − 1
// Always computes for ZINB, Poisson, and Fixed using the same effective mean.

function zinbStats(mu, theta, pi0) {
  const mean     = (1 - pi0) * mu
  const variance = mean * (1 + mu * (pi0 + 1 / theta))
  return { mean, variance }
}

function prestonSibs(mean, variance) {
  return mean > 0 ? mean + variance / mean - 1 : 0
}

export default function PrestonPanel({ params, activeCohort }) {
  const { mu, theta, pi0 } = params

  const zinb    = zinbStats(mu, theta, pi0)
  const poisson = { mean: zinb.mean, variance: zinb.mean }          // same mean, Poisson variance
  const fixed   = { mean: zinb.mean, variance: 0 }                  // same mean, zero variance

  const sibsZINB    = prestonSibs(zinb.mean,    zinb.variance)
  const sibsPoisson = prestonSibs(poisson.mean, poisson.variance)
  const sibsFixed   = prestonSibs(fixed.mean,   fixed.variance)

  const gap = sibsPoisson > 0
    ? ((sibsZINB - sibsPoisson) / sibsPoisson * 100).toFixed(0)
    : '—'

  return (
    <div className="preston-panel">
      <h3 className="preston-title">Preston Decomposition</h3>
      <p className="preston-formula">
        E[siblings] = X̄ + σ²<sub>X</sub>/X̄ − 1
      </p>

      <div className="preston-grid">
        <div className="preston-item">
          <span className="preston-item-label">Fertility mean X̄</span>
          <span className="preston-item-value">{zinb.mean.toFixed(3)}</span>
        </div>
        <div className="preston-item">
          <span className="preston-item-label">ZINB variance σ²</span>
          <span className="preston-item-value">{zinb.variance.toFixed(3)}</span>
        </div>
        <div className="preston-item">
          <span className="preston-item-label">Overdispersion σ²/X̄</span>
          <span className="preston-item-value highlight">{(zinb.variance / zinb.mean).toFixed(3)}</span>
        </div>
        <div className="preston-item">
          <span className="preston-item-label">ZINB → E[siblings]</span>
          <span className="preston-item-value highlight">{sibsZINB.toFixed(3)}</span>
        </div>
      </div>

      <div className="variance-comparison">
        <p className="variance-comparison-title">Expected Siblings by Model</p>
        <table className="variance-table">
          <thead>
            <tr>
              <th>Model</th>
              <th>Variance</th>
              <th>E[siblings]</th>
            </tr>
          </thead>
          <tbody>
            {activeCohort && (
              <tr className="row-empirical">
                <td>Empirical ({activeCohort.year})</td>
                <td>{activeCohort.empVariance.toFixed(3)}</td>
                <td>{activeCohort.empSiblingMean.toFixed(3)}</td>
              </tr>
            )}
            <tr className="row-active">
              <td style={{ color: '#c05621', fontWeight: 600 }}>ZINB</td>
              <td>{zinb.variance.toFixed(3)}</td>
              <td style={{ fontWeight: 700 }}>{sibsZINB.toFixed(3)}</td>
            </tr>
            <tr>
              <td style={{ color: '#3870ae' }}>Poisson</td>
              <td>{poisson.variance.toFixed(3)}</td>
              <td>{sibsPoisson.toFixed(3)}</td>
            </tr>
            <tr>
              <td style={{ color: '#888' }}>Fixed</td>
              <td>0.000</td>
              <td>{sibsFixed.toFixed(3)}</td>
            </tr>
          </tbody>
        </table>
        <div className="fmf-insight" style={{ marginTop: '0.75rem' }}>
          ZINB predicts <strong>{gap}% more siblings</strong> than Poisson — purely
          from overdispersion. This gap compounds at each generational step.
        </div>
      </div>
    </div>
  )
}
