import { useMemo, useRef } from 'react'
import {
  ComposedChart, Bar, Line, XAxis, YAxis,
  CartesianGrid, Tooltip, Legend, ResponsiveContainer,
} from 'recharts'
import {
  siblingPMFFromFertility,
  BUTTERFLY_MAX_K,
} from '../lib/pmfUtils.js'
import { poissonPMF, zinbPMF, nbPMF } from '../lib/distributions.js'
import { relativeError, pearsonChiSq, aic, bestFit } from '../lib/fitMetrics.js'
import { downloadCSV, downloadChartPng, svgFromRef } from '../lib/exportUtils.js'
import GenerationDiagram from './GenerationDiagram.jsx'
import DataProvenance from './DataProvenance.jsx'

const PMF_CSV_COLUMNS = ['count', 'empirical', 'poissonFit', 'zinbFit']
const N_BINS = BUTTERFLY_MAX_K + 1  // 0..11 plus the 12+ bin

function pct(v) { return (Math.abs(v) * 100).toFixed(1) + '%' }

// Show a moment value with its relative error vs empirical, e.g. "3.387 (0.7%)".
function moment(value, err) {
  return err == null ? value.toFixed(3) : `${value.toFixed(3)} (${(err * 100).toFixed(1)}%)`
}

// Format a large fit statistic: small values keep a decimal, large ones get
// thousands separators (χ² and ΔAIC are huge at census sample sizes).
function fmtStat(v) {
  if (v == null) return '—'
  if (v === 0) return '0'
  return v < 100 ? v.toFixed(1) : Math.round(v).toLocaleString()
}

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

function PmfChart({ title, data, xLabel, note, height = 300, csvName }) {
  const ref = useRef(null)
  return (
    <div>
      <div className="fmf-chart-header">
        <h3>{title}</h3>
        <div className="chart-export">
          <button onClick={() => downloadCSV(`${csvName}.csv`, data, PMF_CSV_COLUMNS)}>CSV</button>
          <button onClick={() => downloadChartPng(svgFromRef(ref), `${csvName}.png`)}>PNG</button>
        </div>
      </div>
      <div ref={ref}>
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
      </div>
      {note && <p className="fmf-chart-note">{note}</p>}
    </div>
  )
}

export default function FertilityModelFit({ dataset, selectedYear, onYearChange }) {
  const cohorts = dataset.cohorts
  const cohort  = cohorts.find(c => c.year === selectedYear) ?? cohorts[cohorts.length - 1]
  const fertProbs = dataset.pmfByYear[cohort.year] ?? null

  const {
    empMean, empVariance, mu, theta, pi0,
    empSiblingMean, empSiblingVariance, zinbVariance,
    zinbSiblingMean, zinbSiblingVariance,
  } = cohort
  const effMean = (1 - pi0) * mu

  // ── Fertility chart data ──
  const fertChartData = useMemo(() => {
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
  }, [fertProbs, empMean, mu, theta, pi0])

  // ── Sibling chart data ──
  // Empirical sibling PMF: size-biased from the fertility PMF
  //   P(Y=k) = (k+1)·P(X=k+1) / E[X]
  // ZINB sibling: analytically NB(mu*(θ+1)/θ, θ+1) — π₀ vanishes
  // Poisson sibling: Poisson(empMean) — same distribution as fertility
  const empSibProbs = useMemo(
    () => (fertProbs ? siblingPMFFromFertility(fertProbs, empMean) : null),
    [fertProbs, empMean]
  )

  const sibChartData = useMemo(() => {
    const zinbSibMu    = mu * (theta + 1) / theta
    const zinbSibTheta = theta + 1
    const rows = []
    let pSum = 0, zSum = 0
    for (let k = 0; k < BUTTERFLY_MAX_K; k++) {
      const poisP = poissonPMF(k, empMean)
      const zinbP = nbPMF(k, zinbSibMu, zinbSibTheta)
      pSum += poisP; zSum += zinbP
      rows.push({ count: String(k), empirical: empSibProbs?.[k] ?? null, poissonFit: poisP, zinbFit: zinbP })
    }
    rows.push({
      count:      '12+',
      empirical:  empSibProbs?.[BUTTERFLY_MAX_K] ?? null,
      poissonFit: Math.max(0, 1 - pSum),
      zinbFit:    Math.max(0, 1 - zSum),
    })
    return rows
  }, [empSibProbs, empMean, mu, theta])

  const N = cohort.sampleSize ?? null

  // ── Moment errors (the intuitive discriminator) ──
  const fertMoments = {
    zinb:    { mean: effMean, meanErr: relativeError(effMean, empMean), variance: zinbVariance, varErr: relativeError(zinbVariance, empVariance) },
    poisson: { mean: empMean, meanErr: relativeError(empMean, empMean), variance: empMean,      varErr: relativeError(empMean, empVariance) },
  }
  const sibMoments = {
    zinb:    { mean: zinbSiblingMean, meanErr: relativeError(zinbSiblingMean, empSiblingMean), variance: zinbSiblingVariance, varErr: relativeError(zinbSiblingVariance, empSiblingVariance) },
    poisson: { mean: empMean,         meanErr: relativeError(empMean, empSiblingMean),         variance: empMean,             varErr: relativeError(empMean, empSiblingVariance) },
  }
  // Sibling distribution is a derived transform (no valid sample-based GOF) →
  // judged by the variance error, which is the discriminator.
  const sibBest = bestFit({ zinb: sibMoments.zinb.varErr, poisson: sibMoments.poisson.varErr })

  // ── Whole-distribution goodness of fit on the observed fertility counts ──
  // Pearson χ² + AIC (ΔAIC penalizes ZINB's extra parameters). Only valid for
  // fertility, which is a real sample; requires the cohort's sample size N.
  // Cheap (13-bin loops), so computed each render rather than memoized.
  const fertGOF = (() => {
    if (!N) return { N: null, best: bestFit({ zinb: fertMoments.zinb.varErr, poisson: fertMoments.poisson.varErr }) }
    const zinbProbs = fertChartData.map(r => r.zinbFit)
    const poisProbs = fertChartData.map(r => r.poissonFit)
    const chiZ = pearsonChiSq(fertProbs, zinbProbs, N)
    const chiP = pearsonChiSq(fertProbs, poisProbs, N)
    const aicZ = aic(fertProbs, zinbProbs, N, 3)  // μ, θ, π₀
    const aicP = aic(fertProbs, poisProbs, N, 1)  // λ
    const minAIC = Math.min(aicZ, aicP)
    const dfZ = N_BINS - 1 - 3, dfP = N_BINS - 1 - 1
    return {
      N,
      zinb:    { redChi: chiZ / dfZ, dAIC: aicZ - minAIC },
      poisson: { redChi: chiP / dfP, dAIC: aicP - minAIC },
      best: bestFit({ zinb: aicZ, poisson: aicP }),
    }
  })()

  const varRatio = (empVariance / empMean).toFixed(1)

  return (
    <div className="fmf-container">

      {/* Year selector */}
      <div className="year-pills">
        {cohorts.map(c => (
          <button
            key={c.year}
            className={`year-pill ${c.year === cohort.year ? 'active' : ''}`}
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
            title={`Fertility Distribution — ${cohort.year}`}
            data={fertChartData}
            xLabel="Children ever born"
            csvName={`fertility_fit_${cohort.year}`}
            note={`Grey bars = observed ${dataset.label} data. Orange = ZINB fitted (θ = ${theta.toFixed(2)}, π₀ = ${pi0.toFixed(3)}). Blue dashed = Poisson. Both models match the empirical mean — they differ in variance.`}
            height={300}
          />

          <div className="fmf-divider" />

          {/* Sibling distribution — explicitly the NEXT generation */}
          <div className="fmf-sibling-section">
            <h3 className="fmf-section-heading">Sibling distribution — the next generation</h3>
            <p className="fmf-section-lead">
              The sibling curve does <strong>not</strong> describe the women observed in {cohort.year}.
              It describes the sibship sizes experienced by <strong>their children</strong>, obtained
              by size-biasing the fertility distribution.
            </p>
            <GenerationDiagram year={cohort.year} cohort={cohort.cohort} />
            <PmfChart
              title={`Sibling Distribution — children of the ${cohort.year} cohort`}
              data={sibChartData}
              xLabel="Number of siblings (among the children)"
              csvName={`sibling_fit_${cohort.year}`}
              note={`Grey bars = empirical sibling PMF, size-biased from the fertility data. Orange = ZINB sibling: NB(μ·(θ+1)/θ, θ+1) — π₀ vanishes because childless women contribute no children. Blue dashed = Poisson sibling: Poisson(X̄).`}
              height={300}
            />
          </div>

        </div>

        {/* ── Right: stats + computed fit ── */}
        <div className="fmf-right">

          {/* About this data lives here so provenance sits beside the charts */}
          <DataProvenance dataset={dataset} />

          <div className="fmf-block">
            <h4 className="fmf-block-title">Fertility — model vs data</h4>
            <table className="fmf-table">
              <thead>
                <tr><th>Model</th><th>Mean (err)</th><th>Variance (err)</th></tr>
              </thead>
              <tbody>
                <tr className="fmf-row-emp">
                  <td>Empirical</td>
                  <td>{empMean.toFixed(3)}</td>
                  <td>{empVariance.toFixed(3)}</td>
                </tr>
                <tr className={`fmf-row-zinb ${fertGOF.best === 'zinb' ? 'fit-winner' : ''}`}>
                  <td>ZINB {fertGOF.best === 'zinb' && <span className="fit-good">✓</span>}</td>
                  <td>{moment(fertMoments.zinb.mean, fertMoments.zinb.meanErr)}</td>
                  <td>{moment(fertMoments.zinb.variance, fertMoments.zinb.varErr)}</td>
                </tr>
                <tr className={`fmf-row-pois ${fertGOF.best === 'poisson' ? 'fit-winner' : ''}`}>
                  <td>Poisson {fertGOF.best === 'poisson' && <span className="fit-good">✓</span>}</td>
                  <td>{moment(fertMoments.poisson.mean, fertMoments.poisson.meanErr)}</td>
                  <td>{moment(fertMoments.poisson.variance, fertMoments.poisson.varErr)}</td>
                </tr>
              </tbody>
            </table>
            <p className="fmf-var-note">
              Mean and variance shown with relative error vs the empirical moment.
              Empirical variance is <strong>{varRatio}× the mean</strong> — Poisson forces variance = mean.
            </p>

            {fertGOF.N && (
              <div className="fmf-gof">
                <div className="fmf-gof-title">Whole-distribution fit · N = {fertGOF.N.toLocaleString()}</div>
                <table className="fmf-table fmf-gof-table">
                  <thead>
                    <tr><th></th><th>χ²/df</th><th>ΔAIC</th></tr>
                  </thead>
                  <tbody>
                    <tr className={fertGOF.best === 'zinb' ? 'fit-winner' : ''}>
                      <td>ZINB</td><td>{fmtStat(fertGOF.zinb.redChi)}</td><td>{fmtStat(fertGOF.zinb.dAIC)}</td>
                    </tr>
                    <tr className={fertGOF.best === 'poisson' ? 'fit-winner' : ''}>
                      <td>Poisson</td><td>{fmtStat(fertGOF.poisson.redChi)}</td><td>{fmtStat(fertGOF.poisson.dAIC)}</td>
                    </tr>
                  </tbody>
                </table>
                <p className="fmf-gof-note">
                  Lower = closer. ΔAIC penalizes ZINB's extra parameters; the lower-AIC model is preferred.
                  At census sample sizes both models are formally rejected by χ², so these are comparative, not pass/fail.
                </p>
              </div>
            )}
          </div>

          <div className="fmf-block">
            <h4 className="fmf-block-title">Sibling — model vs data</h4>
            <table className="fmf-table">
              <thead>
                <tr><th>Model</th><th>Mean (err)</th><th>Variance (err)</th></tr>
              </thead>
              <tbody>
                <tr className="fmf-row-emp">
                  <td>Empirical</td>
                  <td>{empSiblingMean.toFixed(3)}</td>
                  <td>{empSiblingVariance.toFixed(3)}</td>
                </tr>
                <tr className={`fmf-row-zinb ${sibBest === 'zinb' ? 'fit-winner' : ''}`}>
                  <td>ZINB → NB {sibBest === 'zinb' && <span className="fit-good">✓</span>}</td>
                  <td>{moment(sibMoments.zinb.mean, sibMoments.zinb.meanErr)}</td>
                  <td>{moment(sibMoments.zinb.variance, sibMoments.zinb.varErr)}</td>
                </tr>
                <tr className={`fmf-row-pois ${sibBest === 'poisson' ? 'fit-winner' : ''}`}>
                  <td>Poisson {sibBest === 'poisson' && <span className="fit-good">✓</span>}</td>
                  <td>{moment(sibMoments.poisson.mean, sibMoments.poisson.meanErr)}</td>
                  <td>{moment(sibMoments.poisson.variance, sibMoments.poisson.varErr)}</td>
                </tr>
              </tbody>
            </table>
            <p className="fmf-var-note">
              The sibling distribution is a size-biased transform of fertility — not an independent
              sample — so it is judged by moment error (closer variance ✓), not a distribution test.
              Under ZINB it is NB(p, s+1): π₀ = {pi0.toFixed(3)} vanishes, since childless women
              contribute no children to the next generation.
            </p>
          </div>

        </div>
      </div>
    </div>
  )
}
