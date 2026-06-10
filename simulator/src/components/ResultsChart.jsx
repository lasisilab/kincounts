import {
  BarChart, Bar, XAxis, YAxis, CartesianGrid,
  Tooltip, Legend, ResponsiveContainer,
} from 'recharts'
import { arrayToPMF, mergePMFs } from '../lib/statsUtils.js'
import { KIN_DISPLAY_CAPS } from '../lib/empiricalData.js'

const KIN_LABELS = {
  children:      'Children',
  siblings:      'Siblings',
  auntsUncles:   'Aunts & Uncles',
  cousins:       'Cousins',
  niecesNephews: 'Nieces & Nephews',
}

const COLORS = {
  children:      '#7bc8a4',
  siblings:      '#4f86c6',
  auntsUncles:   '#e8704a',
  cousins:       '#9b59b6',
  niecesNephews: '#f0c040',
}

function pctFormatter(v) { return (v * 100).toFixed(1) + '%' }

export default function ResultsChart({ results, selectedKin }) {
  if (!results) {
    return (
      <div className="empty-chart">
        <p>Set parameters and click <strong>Run Simulation</strong> to see results.</p>
      </div>
    )
  }

  const active = Object.entries(selectedKin).filter(([, v]) => v).map(([k]) => k)
  if (active.length === 0) return <div className="empty-chart">Select at least one kin type.</div>

  const pmfMap = {}
  for (const key of active) {
    if (results[key]) pmfMap[key] = arrayToPMF(results[key])
  }

  const data = mergePMFs(pmfMap, KIN_DISPLAY_CAPS)

  const capLabels = {}
  for (const key of active) {
    const cap = KIN_DISPLAY_CAPS[key]
    if (cap !== undefined) capLabels[cap] = `${cap}+`
  }
  const xTickFormatter = v => capLabels[v] ?? String(v)

  return (
    <>
      <p className="chart-cap-note">
        X-axis capped at empirical 99.9th-percentile maxima — tail mass folded into final bin.
      </p>
      <ResponsiveContainer width="100%" height={380}>
        <BarChart
          data={data}
          barSize={8}
          barCategoryGap="15%"
          barGap={1}
          margin={{ top: 10, right: 20, left: 10, bottom: 55 }}
        >
          <CartesianGrid strokeDasharray="3 3" vertical={false} stroke="#f0f0f0" />
          <XAxis
            dataKey="count"
            tickFormatter={xTickFormatter}
            label={{ value: 'Number of Relatives', position: 'insideBottom', offset: -38, fontSize: 12 }}
            tick={{ fontSize: 11, angle: -45, textAnchor: 'end', dy: 4 }}
            interval={0}
          />
          <YAxis
            tickFormatter={pctFormatter}
            tick={{ fontSize: 11 }}
            width={50}
            label={{ value: 'Probability', angle: -90, position: 'insideLeft', offset: -5, fontSize: 12 }}
          />
          <Tooltip formatter={(v, name) => [pctFormatter(v), KIN_LABELS[name] ?? name]} />
          <Legend
            verticalAlign="top"
            formatter={key => KIN_LABELS[key] ?? key}
            wrapperStyle={{ fontSize: 12, paddingBottom: 12 }}
          />
          {active.map(key => (
            <Bar key={key} dataKey={key} fill={COLORS[key]} fillOpacity={0.85} isAnimationActive={false} />
          ))}
        </BarChart>
      </ResponsiveContainer>
    </>
  )
}
