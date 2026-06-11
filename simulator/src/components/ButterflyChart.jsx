import {
  BarChart, Bar, XAxis, YAxis, CartesianGrid,
  ReferenceLine, ResponsiveContainer, Tooltip,
} from 'recharts'

const FERT_COLOR = '#4f86c6'
const SIB_COLOR  = '#e8704a'

function absTickFormatter(v) {
  return (Math.abs(v) * 100).toFixed(0) + '%'
}

function CustomTooltip({ active, payload, label }) {
  if (!active || !payload?.length) return null
  const fert = payload.find(p => p.dataKey === 'fert')
  const sib  = payload.find(p => p.dataKey === 'sib')
  return (
    <div className="butterfly-tooltip">
      <p className="bt-label">Count: {label}</p>
      {fert && <p style={{ color: FERT_COLOR }}>Fertility: {(Math.abs(fert.value) * 100).toFixed(1)}%</p>}
      {sib  && <p style={{ color: SIB_COLOR  }}>Sibling:  {(sib.value  * 100).toFixed(1)}%</p>}
    </div>
  )
}

export default function ButterflyChart({ data, year, cohort, isEmpirical }) {
  // Symmetric domain for x-axis
  const maxP = data.reduce((m, d) => Math.max(m, Math.abs(d.fert), d.sib), 0)
  const domainMax = Math.ceil(maxP * 100 + 2) / 100
  const domain = [-domainMax, domainMax]

  return (
    <div className="butterfly-panel">
      <p className="butterfly-year">
        {year}
        {isEmpirical
          ? <span className="bt-badge bt-badge--emp">empirical</span>
          : <span className="bt-badge bt-badge--theory">ZINB</span>
        }
      </p>
      <p className="butterfly-cohort">b. {cohort}</p>
      <ResponsiveContainer width="100%" height={260}>
        <BarChart
          layout="vertical"
          data={data}
          margin={{ top: 4, right: 8, left: 4, bottom: 4 }}
          barSize={9}
          barCategoryGap="20%"
          barGap={0}
        >
          <CartesianGrid horizontal={false} stroke="#f2f2f2" />
          <XAxis
            type="number"
            domain={domain}
            tickFormatter={absTickFormatter}
            tick={{ fontSize: 9 }}
            tickCount={5}
          />
          <YAxis
            type="category"
            dataKey="count"
            tick={{ fontSize: 9 }}
            width={26}
          />
          <ReferenceLine x={0} stroke="#555" strokeWidth={1.5} />
          <Tooltip content={<CustomTooltip />} />
          <Bar dataKey="fert" fill={FERT_COLOR} fillOpacity={0.75} isAnimationActive={false} />
          <Bar dataKey="sib"  fill={SIB_COLOR}  fillOpacity={1.0}  isAnimationActive={false} />
        </BarChart>
      </ResponsiveContainer>
    </div>
  )
}
