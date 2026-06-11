import { IPUMS_COHORTS } from '../lib/empiricalData.js'

function badge(predicted, empirical) {
  if (!empirical) return ''
  const abs = Math.abs((predicted - empirical) / empirical)
  if (abs < 0.02) return 'err-good'
  if (abs < 0.06) return 'err-ok'
  return 'err-bad'
}

function fmt(v) { return v.toFixed(3) }

export default function CohortValidationTable() {
  const rows = IPUMS_COHORTS.map(c => {
    const effectiveMean = (1 - c.pi0) * c.mu  // ZINB effective fertility mean

    // ── Fertility stats ──
    const fert = {
      emp:     { mean: c.empMean,     variance: c.empVariance  },
      zinb:    { mean: effectiveMean, variance: c.zinbVariance },
      poisson: { mean: c.empMean,     variance: c.empMean      },  // σ² = X̄
    }

    // ── Sibling stats ──
    // ZINB sibling stats stored directly from repo's fertility_estimation.csv
    // (avoids rounding errors that arise from analytically recomputing NB(mu*(theta+1)/theta, theta+1)
    //  using the 3-decimal-place rounded mu/theta params stored in IPUMS_COHORTS)
    //
    // Poisson: sibling ~ Poisson(empMean) — same distribution, mean = variance
    const poisSibMu  = c.empMean
    const poisSibVar = c.empMean

    const sib = {
      emp:     { mean: c.empSiblingMean,    variance: c.empSiblingVariance    },
      zinb:    { mean: c.zinbSiblingMean,   variance: c.zinbSiblingVariance   },
      poisson: { mean: poisSibMu,           variance: poisSibVar              },
    }

    return { ...c, fert, sib }
  })

  return (
    <div className="validation-block">
      <div className="validation-header">
        <div>
          <h3 className="validation-title">Fertility &amp; Sibling Distributions — Model vs Empirical</h3>
          <p className="validation-subtitle">
            For each cohort, the ZINB and Poisson models share the same fertility mean.
            Differences in sibling predictions come entirely from how each model handles variance.
            ZINB sibling distribution: NB(μ·(θ+1)/θ, θ+1) — π₀ vanishes from the size-biased transform.
          </p>
        </div>
      </div>

      <table className="validation-table">
        <thead>
          <tr>
            <th rowSpan={2}>Census</th>
            <th rowSpan={2}>Cohort</th>
            <th colSpan={3} className="validation-group-header">Fertility Mean</th>
            <th colSpan={3} className="validation-group-header">Fertility Variance</th>
            <th colSpan={3} className="validation-group-header">Sibling Mean</th>
            <th colSpan={3} className="validation-group-header">Sibling Variance</th>
          </tr>
          <tr>
            <th className="col-emp">Emp</th><th className="col-zinb">ZINB</th><th className="col-pois">Pois</th>
            <th className="col-emp">Emp</th><th className="col-zinb">ZINB</th><th className="col-pois">Pois</th>
            <th className="col-emp">Emp</th><th className="col-zinb">ZINB</th><th className="col-pois">Pois</th>
            <th className="col-emp">Emp</th><th className="col-zinb">ZINB</th><th className="col-pois">Pois</th>
          </tr>
        </thead>
        <tbody>
          {rows.map(r => (
            <tr key={r.year}>
              <td className="col-year">{r.year}</td>
              <td className="col-cohort">{r.cohort}</td>

              {/* Fertility mean */}
              <td className="col-emp">{fmt(r.fert.emp.mean)}</td>
              <td className={`col-zinb ${badge(r.fert.zinb.mean, r.fert.emp.mean)}`}>{fmt(r.fert.zinb.mean)}</td>
              <td className={`col-pois ${badge(r.fert.poisson.mean, r.fert.emp.mean)}`}>{fmt(r.fert.poisson.mean)}</td>

              {/* Fertility variance */}
              <td className="col-emp">{fmt(r.fert.emp.variance)}</td>
              <td className={`col-zinb ${badge(r.fert.zinb.variance, r.fert.emp.variance)}`}>{fmt(r.fert.zinb.variance)}</td>
              <td className={`col-pois ${badge(r.fert.poisson.variance, r.fert.emp.variance)}`}>{fmt(r.fert.poisson.variance)}</td>

              {/* Sibling mean */}
              <td className="col-emp">{fmt(r.sib.emp.mean)}</td>
              <td className={`col-zinb ${badge(r.sib.zinb.mean, r.sib.emp.mean)}`}>{fmt(r.sib.zinb.mean)}</td>
              <td className={`col-pois ${badge(r.sib.poisson.mean, r.sib.emp.mean)}`}>{fmt(r.sib.poisson.mean)}</td>

              {/* Sibling variance */}
              <td className="col-emp">{fmt(r.sib.emp.variance)}</td>
              <td className={`col-zinb ${badge(r.sib.zinb.variance, r.sib.emp.variance)}`}>{fmt(r.sib.zinb.variance)}</td>
              <td className={`col-pois ${badge(r.sib.poisson.variance, r.sib.emp.variance)}`}>{fmt(r.sib.poisson.variance)}</td>
            </tr>
          ))}
        </tbody>
      </table>

      <p className="validation-note">
        Green = within 2% of empirical · Yellow = within 6% · Red = &gt;6% error.
        Poisson fertility variance equals its mean by definition — that single constraint
        propagates into systematic underprediction of both sibling mean and sibling variance.
      </p>
    </div>
  )
}
