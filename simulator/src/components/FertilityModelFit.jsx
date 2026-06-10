import { useMemo } from 'react'
import {
  ComposedChart, Bar, Line, XAxis, YAxis,
  CartesianGrid, Tooltip, Legend, ResponsiveContainer,
} from 'recharts'
import { IPUMS_COHORTS } from '../lib/empiricalData.js'
import {
  parseFertilityCSV,
  siblingPMFFromFertility,
  BUTTERFLY_MAX_K,
} from '../lib/pmfUtils.js'
import { poissonPMF, zinbPMF, nbPMF } from '../lib/distributions.js'

const empiricalRaw = import.meta.glob(
  '/src/data/fertility_pmf_*.csv',
  { eager: true, as: 'raw' }
)

function pct(v) { return (Math.abs(v) * 100).toFixed(1) + '%' }

function CustomTooltip({ active, payload, label }) {
  if (!active || !payload?.length) return null
  return (
    <div className="fmf-tooltip">
      <p className="fmf-tooltip-label">Count: {label}</p>
      {payload.map(p => (
        <p key={p.name} style={{ color: p.color }}>{p.name}: {pct(p.value)}</p>
      ))}
    </div>
  )
}

function PmfChart({ title, badge, data, xLabel, note, height = 300 }) {
  return (
    <div>
      <div className="fmf-chart-header">
        <h3>{title}</h3>
        {badge}
      </div>
      <ResponsiveContainer width="100%" height={height}>
        <ComposedChart data={data} margin={{ top: 8, right: 16, left: 0, bottom: 36 }}>
          <CartesianGrid strokeDasharray="3 3" vertical={false} stroke="#f0f0f0" />
          <XAxis
            dataKey="count"
            tick={{ fontSize: 11 }}
            label={{ value: xLabel, position: 'insideBottom', offset: -20, fontSize: 12 }}
          />
          <YAxis
            tickFormatter={v => (v * 100).toFixed(0) + '%'}
            tick={{ fontSize: 11 }}
            width={44}
            label={{ value: 'Probability', angle: -90, position: 'insideLeft', offset: 8, fontSize: 12 }}
          />
          <Tooltip content={<CustomTooltip />} />
          <Legend verticalAlign="top" wrapperStyle={{ fontSize: 12, paddingBottom: 6 }} />
          <Bar  dataKey="empirical"  name="Empirical"     fill="#94a3b8" fillOpacity={0.55} isAnimationActive={false} />
          <Line dataKey="zinbFit"    name="ZINB (fitted)"  stroke="#e8704a" strokeWidth={2.5} dot={{ r: 3 }} type="monotone" isAnimationActive={false} />
          <Line dataKey="poissonFit" name="Poisson"        stroke="#4f86c6" strokeWidth={2} strokeDasharray="6 3" dot={{ r: 3 }} type="monotone" isAnimationActive={false} />
        </ComposedChart>
      </ResponsiveContainer>
      {note && <p className="fmf-chart-note">{note}</p>}
    </div>
  )
}

export default function FertilityModelFit({ selectedYear, onYearChange }) {
  const cohort  = IPUMS_COHORTS.find(c => c.year === selectedYear) ?? IPUMS_COHORTS[4]
  const rawCSV  = empiricalRaw[`/src/data/fertility_pmf_${selectedYear}.csv`]

  const fertProbs = useMemo(() =>
    rawCSV ? parseFertilityCSV(rawCSV) : null
  , [rawCSV])

  // ── Fertility chart data ──
  const fertChartData = useMemo(() => {
    const { empMean, mu, theta, pi0 } = cohort
    const rows = []
    let pSum = 0, zSum = 0
    for (let k = 0; k < BUTTERFLY_MAX_K; k++) {
      const poisP = poissonPMF(k, empMean)
      const zinbP = zinbPMF(k, mu, theta, pi0)
      pSum += poisP; zSum += zinbP
      rows.push({ count: String(k), empirical: fertProbs?.[k] ?? null, poissonFit: poisP, zinbFit: zinbP })
    }
    rows.push({
      count: '12+',
      empirical:  fertProbs?.[BUTTERFLY_MAX_K] ?? null,
      poissonFit: Math.max(0, 1 - pSum),
      zinbFit:    Math.max(0, 1 - zSum),
    })
    return rows
  }, [fertProbs, cohort])

  // ── Sibling chart data ──
  // Empirical sibling PMF: size-biased from fertility CSV
  //   P(Y=k) = (k+1)·P(X=k+1) / E[X]
  // ZINB sibling: analytically NB(mu*(θ+1)/θ, θ+1) — π₀ vanishes
  // Poisson sibling: Poisson(empMean) — same distribution as fertility
  const sibChartData = useMemo(() => {
    const { empMean, mu, theta } = cohort
    const empSibProbs = fertProbs ? siblingPMFFromFertility(fertProbs, empMean) : null

    const zinbSibMu    = mu * (theta + 1) / theta
    const zinbSibTheta = theta + 1

    const rows = []
    let pSum = 0, zSum = 0
    for (let k = 0; k < BUTTERFLY_MAX_K; k++) {
      const poisP = poissonPMF(k, empMean)
      const zinbP = nbPMF(k, zinbSibMu, zinbSibTheta)
      pSum += poisP; zSum += zinbP
      rows.push({
        count:      String(k),
        empirical:  empSibProbs?.[k] ?? null,
        poissonFit: poisP,
        zinbFit:    zinbP,
      })
    }
    rows.push({
      count:      '12+',
      empirical:  empSibProbs?.[BUTTERFLY_MAX_K] ?? null,
      poissonFit: Math.max(0, 1 - pSum),
      zinbFit:    Math.max(0, 1 - zSum),
    })
    return rows
  }, [fertProbs, cohort])

  // ── Stats ──
  // zinbSiblingMean/Variance stored directly from repo CSV (avoids rounding errors
  // from recomputing NB(mu*(θ+1)/θ, θ+1) with 3-decimal-place rounded mu/theta)
  const {
    empMean, empVariance, empSiblingMean, empSiblingVariance,
    zinbVariance, zinbSiblingMean, zinbSiblingVariance,
  } = cohort
  const varRatio = (empVariance / empMean).toFixed(1)

  return (
    <div className="fmf-container">

      {/* Year selector */}
      <div className="year-pills">
        {IPUMS_COHORTS.map(c => (
          <button
            key={c.year}
            className={`year-pill ${c.year === selectedYear ? 'active' : ''}`}
            onClick={() => onYearChange(c.year)}
          >
            <span className="year-pill-year">{c.year}</span>
            <span className="year-pill-cohort">{c.cohort}</span>
          </button>
        ))}
      </div>

      <div className="fmf-body">

        {/* ── Left: fertility chart + sibling chart ── */}
        <div className="fmf-chart-wrap">

          <PmfChart
            title={`Fertility Distribution — ${selectedYear}`}
            data={fertChartData}
            xLabel="Children ever born"
            note={`Grey bars = observed IPUMS census data. Orange = ZINB fitted (θ = ${cohort.theta.toFixed(2)}, π₀ = ${cohort.pi0.toFixed(3)}). Blue dashed = Poisson. Both models share the same mean — the difference is entirely in variance.`}
            height={300}
          />

          <div className="fmf-divider" />

          <PmfChart
            title={`Sibling Distribution — ${selectedYear}`}
            data={sibChartData}
            xLabel="Number of siblings"
            note={`Grey bars = empirical sibling PMF derived via P(Y=k) = (k+1)·P(X=k+1)/E[X] applied to IPUMS data. Orange = ZINB sibling: NB(p, s+1) — π₀ vanishes from the sibling distribution. Blue dashed = Poisson sibling: Poisson(X̄).`}
            height={300}
          />

        </div>

        {/* ── Right: stats panels ── */}
        <div className="fmf-right">

          <div className="fmf-block">
            <h4 className="fmf-block-title">Fertility</h4>
            <table className="fmf-table">
              <thead>
                <tr><th>Model</th><th>Mean</th><th>Variance</th><th></th></tr>
              </thead>
              <tbody>
                <tr className="fmf-row-emp">
                  <td>Empirical</td>
                  <td>{empMean.toFixed(3)}</td>
                  <td>{empVariance.toFixed(3)}</td>
                  <td>—</td>
                </tr>
                <tr className="fmf-row-zinb">
                  <td>ZINB</td>
                  <td>{empMean.toFixed(3)}</td>
                  <td>{zinbVariance.toFixed(3)}</td>
                  <td className="fit-good">✓</td>
                </tr>
                <tr className="fmf-row-pois">
                  <td>Poisson</td>
                  <td>{empMean.toFixed(3)}</td>
                  <td>{empMean.toFixed(3)}</td>
                  <td className="fit-bad">✗</td>
                </tr>
              </tbody>
            </table>
            <p className="fmf-var-note">
              Empirical variance is <strong>{varRatio}× larger</strong> than the mean.
              Poisson cannot capture this overdispersion.
            </p>
          </div>

          <div className="fmf-block">
            <h4 className="fmf-block-title">Sibling</h4>
            <table className="fmf-table">
              <thead>
                <tr><th>Model</th><th>Mean</th><th>Variance</th><th></th></tr>
              </thead>
              <tbody>
                <tr className="fmf-row-emp">
                  <td>Empirical</td>
                  <td>{empSiblingMean.toFixed(3)}</td>
                  <td>{empSiblingVariance.toFixed(3)}</td>
                  <td>—</td>
                </tr>
                <tr className="fmf-row-zinb">
                  <td>ZINB → NB</td>
                  <td>{zinbSiblingMean.toFixed(3)}</td>
                  <td>{zinbSiblingVariance.toFixed(3)}</td>
                  <td className="fit-good">✓</td>
                </tr>
                <tr className="fmf-row-pois">
                  <td>Poisson</td>
                  <td>{empMean.toFixed(3)}</td>
                  <td>{empMean.toFixed(3)}</td>
                  <td className="fit-bad">✗</td>
                </tr>
              </tbody>
            </table>
            <p className="fmf-var-note">
              Under ZINB fertility, the sibling distribution is NB(p, s+1).
              Zero-inflation π₀ = {cohort.pi0.toFixed(3)} disappears — childless
              women contribute no individuals to the next generation.
            </p>
          </div>

        </div>
      </div>
    </div>
  )
}
