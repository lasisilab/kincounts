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
import { poissonPMF, nbPMF } from '../lib/distributions.js'

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
        <p key={p.name} style={{ color: p.color ?? p.fill }}>
          {p.name}: {pct(p.value)}
        </p>
      ))}
    </div>
  )
}

function PmfChart({ title, data, badge }) {
  return (
    <div className="sbp-chart-wrap">
      <div className="fmf-chart-header">
        <h4 className="sbp-chart-title">{title}</h4>
        {badge}
      </div>
      <ResponsiveContainer width="100%" height={300}>
        <ComposedChart data={data} margin={{ top: 6, right: 12, left: 0, bottom: 36 }}>
          <CartesianGrid strokeDasharray="3 3" vertical={false} stroke="#f0f0f0" />
          <XAxis
            dataKey="count"
            tick={{ fontSize: 10 }}
            label={{ value: 'Count', position: 'insideBottom', offset: -22, fontSize: 11 }}
          />
          <YAxis
            tickFormatter={v => (v * 100).toFixed(0) + '%'}
            tick={{ fontSize: 10 }}
            width={40}
          />
          <Tooltip content={<CustomTooltip />} />
          <Legend verticalAlign="top" wrapperStyle={{ fontSize: 11, paddingBottom: 4 }} />
          <Bar  dataKey="empirical" name="Empirical"    fill="#94a3b8" fillOpacity={0.55} isAnimationActive={false} />
          <Line dataKey="zinb"      name="ZINB"         stroke="#e8704a" strokeWidth={2.5} dot={{ r: 2.5 }} type="monotone" isAnimationActive={false} />
          <Line dataKey="poisson"   name="Poisson"      stroke="#4f86c6" strokeWidth={2} strokeDasharray="6 3" dot={{ r: 2.5 }} type="monotone" isAnimationActive={false} />
        </ComposedChart>
      </ResponsiveContainer>
    </div>
  )
}

export default function SizeBiasingPanel({ selectedYear }) {
  const cohort = IPUMS_COHORTS.find(c => c.year === selectedYear) ?? IPUMS_COHORTS[4]
  const { mu, theta, pi0, empMean } = cohort

  const rawCSV = empiricalRaw[`/src/data/fertility_pmf_${selectedYear}.csv`]

  const { fertData, sibData } = useMemo(() => {
    // Empirical fertility PMF
    const empFertProbs = rawCSV ? parseFertilityCSV(rawCSV) : null

    // ZINB sibling distribution is analytically NB(mu*(theta+1)/theta, theta+1)
    // — π₀ disappears entirely from the sibling distribution
    const zinbSibMu    = mu * (theta + 1) / theta
    const zinbSibTheta = theta + 1

    const fertRows = []
    const sibRows  = []
    let pFertSum = 0, zFertSum = 0
    let pSibSum  = 0, zSibSum  = 0

    for (let k = 0; k < BUTTERFLY_MAX_K; k++) {
      const poissonFert = poissonPMF(k, empMean)
      const zinbFert    = empFertProbs ? null : null  // bars = empirical only
      pFertSum += poissonFert

      const poissonSib = poissonPMF(k, empMean)  // Poisson sibling = Poisson(empMean)
      const zinbSib    = nbPMF(k, zinbSibMu, zinbSibTheta)
      pSibSum  += poissonSib
      zSibSum  += zinbSib

      fertRows.push({
        count:    String(k),
        empirical: empFertProbs?.[k] ?? null,
        zinb:     null,   // fitted line not shown on fertility side (already in Tab 1 chart)
        poisson:  poissonFert,
      })
      sibRows.push({
        count:    String(k),
        empirical: empFertProbs ? siblingPMFFromFertility(empFertProbs, empMean)[k] : null,
        zinb:     zinbSib,
        poisson:  poissonSib,
      })
    }

    // 12+ tail bins
    fertRows.push({
      count:    '12+',
      empirical: empFertProbs ? empFertProbs[BUTTERFLY_MAX_K] : null,
      zinb:     null,
      poisson:  Math.max(0, 1 - pFertSum),
    })
    sibRows.push({
      count:    '12+',
      empirical: empFertProbs
        ? siblingPMFFromFertility(empFertProbs, empMean)[BUTTERFLY_MAX_K]
        : null,
      zinb:    Math.max(0, 1 - zSibSum),
      poisson: Math.max(0, 1 - pSibSum),
    })

    return { fertData: fertRows, sibData: sibRows }
  }, [rawCSV, cohort])

  const zinbSibMean = (mu * (theta + 1) / theta).toFixed(3)

  return (
    <div className="sbp-container">
      <div className="sbp-header">
        <h3 className="sbp-title">Size-Biasing Transformation — {selectedYear}</h3>
        <p className="sbp-subtitle">
          Applying P(Y = k) = (k + 1) · P(X = k + 1) / E[X] to the fertility distribution
          yields the sibling distribution. The charts below show this shift for the {selectedYear} cohort.
        </p>
      </div>

      <div className="sbp-charts">
        <PmfChart
          title="Fertility Distribution"
          data={fertData}
          badge={<span className="fmf-badge fmf-badge--emp">input</span>}
        />
        <div className="sbp-arrow">→</div>
        <PmfChart
          title="Sibling Distribution (derived)"
          data={sibData}
          badge={<span className="fmf-badge fmf-badge--theory">size-biased</span>}
        />
      </div>

      {/* Key analytical result callout */}
      <div className="sbp-callout">
        <div className="sbp-callout-formula">
          ZINB(p, s, π₀) → <strong>NB(p, s + 1)</strong>
        </div>
        <p className="sbp-callout-text">
          Under ZINB fertility, the induced sibling distribution is a pure Negative Binomial —
          the zero-inflation probability <strong>π₀ disappears entirely</strong>.
          Women with no children contribute no individuals to the child generation,
          so their childlessness is invisible to siblings.
          For the {selectedYear} cohort, the ZINB sibling mean is <strong>{zinbSibMean}</strong>,
          compared to the Poisson prediction of <strong>{Number(cohort.empMean).toFixed(3)}</strong>.
        </p>
      </div>
    </div>
  )
}
