import { useMemo } from 'react'
import ButterflyChart from './ButterflyChart.jsx'
import { IPUMS_COHORTS } from '../lib/empiricalData.js'
import {
  fertilityPMFArray,
  siblingPMFFromFertility,
  parseFertilityCSV,
  buildButterflyData,
} from '../lib/pmfUtils.js'

// Attempt to load empirical CSV files placed in src/data/
// Files must be named: fertility_pmf_1950.csv … fertility_pmf_1990.csv
// If the files aren't present, the glob returns an empty object and the
// component falls back to the ZINB theoretical PMF automatically.
const empiricalRaw = import.meta.glob(
  '/src/data/fertility_pmf_*.csv',
  { eager: true, query: '?raw', import: 'default' }
)

function keyForYear(year) {
  return `/src/data/fertility_pmf_${year}.csv`
}

export default function ButterflyPanel() {
  const panels = useMemo(() => {
    return IPUMS_COHORTS.map(cohort => {
      const rawCSV = empiricalRaw[keyForYear(cohort.year)]
      let fertProbs, isEmpirical

      if (rawCSV) {
        fertProbs  = parseFertilityCSV(rawCSV)
        isEmpirical = true
      } else {
        fertProbs  = fertilityPMFArray(cohort.mu, cohort.theta, cohort.pi0)
        isEmpirical = false
      }

      const fertMean = isEmpirical ? cohort.empMean : cohort.empMean
      const sibProbs = siblingPMFFromFertility(fertProbs, fertMean)
      const data     = buildButterflyData(fertProbs, sibProbs)

      return { cohort, fertProbs, sibProbs, data, isEmpirical }
    })
  }, [])  // recomputes only on mount; CSVs are bundled at build time

  const anyEmpirical = panels.some(p => p.isEmpirical)

  return (
    <div className="butterfly-section">
      <h2>Fertility ↔ Sibling PMFs by Cohort</h2>
      <div className="butterfly-legend">
        <span className="bt-legend-item bt-legend--fert">◀ Fertility</span>
        <span className="bt-legend-item bt-legend--sib">Sibling ▶</span>
        {!anyEmpirical && (
          <span className="bt-legend-note">
            Showing ZINB theoretical — drop <code>fertility_pmf_YYYY.csv</code> files into <code>src/data/</code> for empirical
          </span>
        )}
      </div>
      <div className="butterfly-grid">
        {panels.map(({ cohort, data, isEmpirical }) => (
          <ButterflyChart
            key={cohort.year}
            data={data}
            year={cohort.year}
            cohort={cohort.cohort}
            isEmpirical={isEmpirical}
          />
        ))}
      </div>
    </div>
  )
}
