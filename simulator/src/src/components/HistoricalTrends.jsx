import {
  LineChart, Line, XAxis, YAxis, CartesianGrid,
  Tooltip, Legend, ResponsiveContainer, ReferenceDot,
} from 'recharts'
import { IPUMS_COHORTS } from '../lib/empiricalData.js'

const data = IPUMS_COHORTS.map(c => ({
  year: c.year,
  'Fertility mean': c.empMean,
  'Sibling mean': c.empSiblingMean,
}))

function CustomTooltip({ active, payload, label }) {
  if (!active || !payload?.length) return null
  const cohort = IPUMS_COHORTS.find(c => c.year === label)
  return (
    <div className="trends-tooltip">
      <p className="tooltip-year">{label} (born {cohort?.cohort})</p>
      {payload.map(p => (
        <p key={p.name} style={{ color: p.color }}>
          {p.name}: {p.value.toFixed(3)}
        </p>
      ))}
      {payload.length === 2 && (
        <p className="tooltip-gap">
          Gap (σ²/X̄ − 1): {(payload[1].value - payload[0].value).toFixed(3)}
        </p>
      )}
    </div>
  )
}

export default function HistoricalTrends({ activeCohort }) {
  return (
    <div className="trends-panel">
      <h2>Historical Cohort Trends (IPUMS Census)</h2>
      <p className="trends-subtitle">
        The sibling mean always exceeds the fertility mean by the overdispersion term σ²<sub>X</sub>/X̄ − 1.
        The gap narrows as cohorts become less overdispersed.
      </p>
      <ResponsiveContainer width="100%" height={280}>
        <LineChart
          data={data}
          margin={{ top: 10, right: 30, left: 0, bottom: 10 }}
        >
          <CartesianGrid strokeDasharray="3 3" stroke="#f0f0f0" />
          <XAxis
            dataKey="year"
            tick={{ fontSize: 12 }}
            label={{ value: 'Census year', position: 'insideBottom', offset: -4, fontSize: 12 }}
          />
          <YAxis
            domain={[1.5, 5]}
            tick={{ fontSize: 12 }}
            label={{ value: 'Mean count', angle: -90, position: 'insideLeft', offset: 10, fontSize: 12 }}
          />
          <Tooltip content={<CustomTooltip />} />
          <Legend
            verticalAlign="top"
            wrapperStyle={{ fontSize: 12, paddingBottom: 6 }}
          />
          <Line
            type="monotone"
            dataKey="Fertility mean"
            stroke="#4f86c6"
            strokeWidth={2.5}
            dot={{ r: 4 }}
            activeDot={{ r: 6 }}
          />
          <Line
            type="monotone"
            dataKey="Sibling mean"
            stroke="#e8704a"
            strokeWidth={2.5}
            dot={{ r: 4 }}
            activeDot={{ r: 6 }}
          />
          {activeCohort && (
            <>
              <ReferenceDot
                x={activeCohort.year}
                y={activeCohort.empMean}
                r={7}
                fill="#4f86c6"
                stroke="#fff"
                strokeWidth={2}
              />
              <ReferenceDot
                x={activeCohort.year}
                y={activeCohort.empSiblingMean}
                r={7}
                fill="#e8704a"
                stroke="#fff"
                strokeWidth={2}
              />
            </>
          )}
        </LineChart>
      </ResponsiveContainer>
    </div>
  )
}
