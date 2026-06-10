import { computeStats } from '../lib/statsUtils.js'

const KIN_LABELS = {
  children:      'Children',
  siblings:      'Siblings',
  auntsUncles:   'Aunts & Uncles',
  cousins:       'Cousins',
  niecesNephews: 'Nieces & Nephews',
}

function fmt(v) { return v.toFixed(3) }
function pct(v) { return (v * 100).toFixed(1) + '%' }

export default function SummaryStats({ results, selectedKin }) {
  if (!results) return null

  const active = Object.entries(selectedKin).filter(([, v]) => v).map(([k]) => k)
  if (active.length === 0) return null

  return (
    <table className="stats-table">
      <thead>
        <tr>
          <th>Kin Type</th>
          <th>Mean</th>
          <th>Variance</th>
          <th>P(0)</th>
        </tr>
      </thead>
      <tbody>
        {active.map(key => {
          const arr = results[key]
          if (!arr) return null
          const { mean, variance, pZero } = computeStats(arr)
          return (
            <tr key={key}>
              <td>{KIN_LABELS[key]}</td>
              <td>{fmt(mean)}</td>
              <td>{fmt(variance)}</td>
              <td>{pct(pZero)}</td>
            </tr>
          )
        })}
      </tbody>
    </table>
  )
}
