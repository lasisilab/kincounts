import {
  LineChart, Line, XAxis, YAxis, CartesianGrid,
  Tooltip, Legend, ResponsiveContainer, ReferenceDot,
} from 'recharts'

const FERT_KEY = 'Completed fertility (mothers)'
const SIB_KEY  = 'Sibship size (their children)'

function CustomTooltip({ active, payload, label, cohortByYear }) {
  if (!active || !payload?.length) return null
  const cohort = cohortByYear[label]
  return (
    <div className="trends-tooltip">
      <p className="tooltip-year">
        {label}{cohort?.cohort ? ` · mothers born ${cohort.cohort}` : ''}
      </p>
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

// Two-line x-axis tick: census year on top, the mothers' birth cohort below it,
// so the chart shows at a glance which birth cohort each census represents.
function CohortTick({ x, y, payload, cohortByYear }) {
  const cohort = cohortByYear[payload.value]
  return (
    <g transform={`translate(${x},${y})`}>
      <text x={0} y={0} dy={14} textAnchor="middle" fontSize={12} fill="#333">
        {payload.value}
      </text>
      {cohort?.cohort && (
        <text x={0} y={0} dy={29} textAnchor="middle" fontSize={10} fill="#888">
          b. {cohort.cohort}
        </text>
      )}
    </g>
  )
}

export default function HistoricalTrends({ dataset, activeCohort }) {
  const cohorts = [...dataset.cohorts].sort((a, b) => a.year - b.year)
  const cohortByYear = Object.fromEntries(cohorts.map(c => [c.year, c]))
  const data = cohorts.map(c => ({
    year: c.year,
    [FERT_KEY]: c.empMean,
    [SIB_KEY]:  c.empSiblingMean,
  }))

  return (
    <div className="trends-panel">
      <h2>Historical Cohort Trends — {dataset.label}</h2>
      <p className="trends-subtitle">
        Each census year is plotted against the mothers' birth cohort shown below it.
        The two lines describe <strong>different generations</strong>: completed fertility is
        measured on the women observed at that census; sibship size is what their children
        experienced. The sibling line always sits above the fertility line by the
        overdispersion term σ²<sub>X</sub>/X̄ − 1, and the gap narrows as cohorts become less overdispersed.
      </p>
      <ResponsiveContainer width="100%" height={300}>
        <LineChart
          data={data}
          margin={{ top: 10, right: 30, left: 0, bottom: 28 }}
        >
          <CartesianGrid strokeDasharray="3 3" stroke="#f0f0f0" />
          <XAxis
            dataKey="year"
            height={48}
            tick={<CohortTick cohortByYear={cohortByYear} />}
            label={{ value: 'Census year (mothers’ birth cohort below)', position: 'insideBottom', offset: -2, fontSize: 12 }}
          />
          <YAxis
            domain={[1.5, 5]}
            tick={{ fontSize: 12 }}
            label={{ value: 'Mean count', angle: -90, position: 'insideLeft', offset: 10, fontSize: 12 }}
          />
          <Tooltip content={<CustomTooltip cohortByYear={cohortByYear} />} />
          <Legend
            verticalAlign="top"
            wrapperStyle={{ fontSize: 12, paddingBottom: 6 }}
          />
          <Line
            type="monotone"
            dataKey={FERT_KEY}
            stroke="#4f86c6"
            strokeWidth={2.5}
            dot={{ r: 4 }}
            activeDot={{ r: 6 }}
          />
          <Line
            type="monotone"
            dataKey={SIB_KEY}
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
