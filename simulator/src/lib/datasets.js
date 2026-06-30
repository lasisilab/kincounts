// Dataset abstraction for the Kincounts simulator.
//
// A "dataset" bundles everything the Fertility Fit view and the Simulator need:
// per-cohort moments + fitted parameters, the empirical fertility PMF per year,
// and provenance text. The app holds one ACTIVE dataset in state. IPUMS is the
// default; an imported dataset has the same shape and fully REPLACES IPUMS, so
// every view and the simulator operate on the user's data with no special-casing.
//
// Cohort object shape (matches empiricalData.js IPUMS_COHORTS):
//   { year, cohort, empMean, empVariance, mu, theta, pi0,
//     empSiblingMean, empSiblingVariance, zinbVariance,
//     zinbSiblingMean, zinbSiblingVariance }
//
// Dataset shape:
//   { id, label, provenance: {...}, cohorts: [...], pmfByYear: { year: probs[] } }

import { IPUMS_COHORTS } from './empiricalData.js'
import { parseFertilityCSV } from './pmfUtils.js'

// Build the IPUMS empirical PMF lookup from the bundled CSVs at build time.
const ipumsRaw = import.meta.glob(
  '/src/data/fertility_pmf_*.csv',
  { eager: true, query: '?raw', import: 'default' }
)
const ipumsPmfByYear = {}
for (const [path, raw] of Object.entries(ipumsRaw)) {
  const m = path.match(/fertility_pmf_(\d+)\.csv$/)
  if (m) ipumsPmfByYear[Number(m[1])] = parseFertilityCSV(raw)
}

export const IPUMS_DATASET = {
  id: 'ipums',
  label: 'IPUMS USA Census',
  provenance: {
    source: 'IPUMS USA (Ruggles et al.), decennial census microdata',
    variable: 'CHBORN — children ever born (excludes stillbirths, adopted, and stepchildren)',
    sample: 'Women aged 50–59 at each census, 1950–1990 — chosen to capture completed fertility while limiting mortality selection',
    cohorts: 'Decade birth cohorts, 1891–1940',
    caveat: 'In the earliest censuses the children-ever-born question covered only ever-married women (1950–1960). Early-cohort fertility therefore excludes never-married mothers; see the methods page.',
    // Rendered Quarto methods page (GitHub Pages). Update if the Pages URL differs.
    methodsHref: 'https://lasisilab.github.io/kincounts/analysis/preprocess.html',
    methodsLabel: 'Preprocessing IPUMS fertility data (methods)',
  },
  cohorts: IPUMS_COHORTS,
  pmfByYear: ipumsPmfByYear,
}

// Per-model default parameters for one generation, derived from a cohort.
// All three models are seeded from the same cohort so switching models keeps the
// generation consistent. Nothing is hard-coded: values come from the dataset.
export function cohortToGenParams(cohort) {
  return {
    zinb:    { mu: cohort.mu, theta: cohort.theta, pi0: cohort.pi0 },
    poisson: { mu: cohort.empMean },
    fixed:   { fertMean: cohort.empMean, sibMean: cohort.empSiblingMean },
  }
}

// Map a dataset's cohorts onto the three simulator generations:
//   focal = newest cohort, grandparent = oldest, parent = the middle cohort.
// Returns { zinb, poisson, fixed } each { focal, parent, grandparent }, plus the
// seeding cohorts so the UI can show which cohort fed each generation.
export function defaultsFromDataset(dataset) {
  const cohorts = [...dataset.cohorts].sort((a, b) => a.year - b.year)
  const n = cohorts.length
  const oldest = cohorts[0]
  const newest = cohorts[n - 1]
  const middle = cohorts[Math.floor((n - 1) / 2)]

  const f = cohortToGenParams(newest)
  const p = cohortToGenParams(middle)
  const g = cohortToGenParams(oldest)

  return {
    byModel: {
      zinb:    { focal: f.zinb,    parent: p.zinb,    grandparent: g.zinb },
      poisson: { focal: f.poisson, parent: p.poisson, grandparent: g.poisson },
      fixed:   { focal: f.fixed,   parent: p.fixed,   grandparent: g.fixed },
    },
    seeds: { focal: newest, parent: middle, grandparent: oldest },
  }
}
