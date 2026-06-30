// Parse user-supplied fertility data (raw counts) and build a dataset object
// with the same shape as IPUMS_DATASET, so imported data REPLACES IPUMS across
// the whole app. The app fits Poisson + ZINB in-browser (see fitZinb.js); the
// user only supplies the empirical distribution.
//
// Expected CSV columns (header row, case-insensitive):
//   year      — census/observation year (integer)
//   label     — birth-cohort label for display (optional, e.g. "1931–1940")
//   children  — number of children ever born (0,1,2,…; the top value is "or more")
//   count     — number of women with that many children
//
// One row per (year, children). Multiple years are tracked as separate cohorts.

import { siblingPMFFromFertility, BUTTERFLY_MAX_K } from './pmfUtils.js'
import { fitZINB } from './fitZinb.js'

const TAIL = BUTTERFLY_MAX_K // 12

function findCol(header, names) {
  for (const n of names) {
    const i = header.indexOf(n)
    if (i !== -1) return i
  }
  return -1
}

export function parseCountsCSV(text) {
  const lines = text.trim().split(/\r?\n/).filter(l => l.trim() !== '')
  if (lines.length < 2) return { error: 'Need a header row and at least one data row.' }

  const header = lines[0].split(',').map(h => h.trim().toLowerCase())
  const iYear  = findCol(header, ['year'])
  const iLabel = findCol(header, ['label', 'cohort'])
  const iChild = findCol(header, ['children', 'child', 'k', 'parity'])
  const iCount = findCol(header, ['count', 'n', 'women', 'frequency', 'freq'])

  if (iYear === -1)  return { error: 'Missing required "year" column.' }
  if (iChild === -1) return { error: 'Missing required "children" column.' }
  if (iCount === -1) return { error: 'Missing required "count" column.' }

  const rows = []
  for (let r = 1; r < lines.length; r++) {
    const cells = lines[r].split(',')
    const year = Number(cells[iYear])
    const children = Number(cells[iChild])
    const count = Number(cells[iCount])
    if (!Number.isFinite(year) || !Number.isFinite(children) || !Number.isFinite(count)) {
      return { error: `Row ${r + 1}: year, children, and count must all be numbers.` }
    }
    if (children < 0 || count < 0) {
      return { error: `Row ${r + 1}: children and count must be non-negative.` }
    }
    const label = iLabel !== -1 ? (cells[iLabel] ?? '').trim() : ''
    rows.push({ year, label: label || null, children, count })
  }
  return { rows }
}

// Aggregate parsed rows into per-year cohorts with exact moments + binned PMF.
function aggregateByYear(rows) {
  const byYear = new Map()
  for (const { year, label, children, count } of rows) {
    if (!byYear.has(year)) byYear.set(year, { year, label: null, pairs: [] })
    const g = byYear.get(year)
    if (label && !g.label) g.label = label
    g.pairs.push([children, count])
  }

  const cohorts = []
  for (const g of byYear.values()) {
    let N = 0, sum = 0
    for (const [k, c] of g.pairs) { N += c; sum += k * c }
    if (N === 0) continue
    const mean = sum / N
    let varSum = 0
    for (const [k, c] of g.pairs) varSum += c * (k - mean) ** 2
    const variance = varSum / N

    const counts = new Array(TAIL + 1).fill(0)
    for (const [k, c] of g.pairs) {
      const bin = k >= TAIL ? TAIL : k
      counts[bin] += c
    }
    const probs = counts.map(c => c / N)

    cohorts.push({
      year: g.year,
      label: g.label || String(g.year),
      N, mean, variance, counts, probs,
    })
  }
  cohorts.sort((a, b) => a.year - b.year)
  return cohorts
}

function pmfMoments(probs) {
  let mean = 0
  for (let k = 0; k < probs.length; k++) mean += k * probs[k]
  let v = 0
  for (let k = 0; k < probs.length; k++) v += probs[k] * (k - mean) ** 2
  return { mean, variance: v }
}

/**
 * Build a dataset object from a raw-counts CSV string.
 * Returns { dataset, warnings } or { error }.
 */
export function buildUserDataset(text, name = 'Your data') {
  const parsed = parseCountsCSV(text)
  if (parsed.error) return { error: parsed.error }

  const agg = aggregateByYear(parsed.rows)
  if (agg.length === 0) return { error: 'No usable rows found (all counts were zero?).' }

  const warnings = []
  const pmfByYear = {}
  const cohorts = agg.map(c => {
    if (c.N < 100) warnings.push(`Year ${c.year}: only ${c.N} women — fits may be unstable.`)

    const { mu, theta, pi0 } = fitZINB(c.counts, c.mean, c.variance)

    // Empirical sibling distribution (size-biased from the fertility PMF).
    const sibProbs = siblingPMFFromFertility(c.probs, c.mean)
    const sib = pmfMoments(sibProbs)

    // ZINB-implied moments (analytic).
    const zinbVariance = (1 - pi0) * (mu + (mu * mu) / theta + pi0 * mu * mu)
    const zinbSibMean = mu * (theta + 1) / theta
    const zinbSibVar  = zinbSibMean + (zinbSibMean * zinbSibMean) / (theta + 1)

    pmfByYear[c.year] = c.probs
    return {
      year: c.year,
      cohort: c.label,
      empMean: c.mean,
      empVariance: c.variance,
      mu, theta, pi0,
      empSiblingMean: sib.mean,
      empSiblingVariance: sib.variance,
      zinbVariance,
      zinbSiblingMean: zinbSibMean,
      zinbSiblingVariance: zinbSibVar,
      sampleSize: c.N,
    }
  })

  const dataset = {
    id: 'user',
    label: name,
    provenance: {
      source: 'User-provided data (imported into the Kincounts simulator)',
      variable: 'Children ever born, as supplied',
      sample: `${cohorts.length} cohort${cohorts.length === 1 ? '' : 's'}; ${cohorts.map(c => c.year).join(', ')}`,
      caveat: 'Poisson and ZINB were fitted in-browser to your binned counts. Compare each model’s discrepancy in the Fertility Fit table to judge fit quality.',
    },
    cohorts,
    pmfByYear,
  }
  return { dataset, warnings }
}

// Produce a fillable raw-counts template from an existing dataset's PMFs,
// scaled to a nominal sample size so users see a real, valid example.
export function makeTemplateCSV(dataset, nominalN = 10000) {
  const lines = ['year,label,children,count']
  const cohorts = [...dataset.cohorts].sort((a, b) => a.year - b.year)
  for (const c of cohorts) {
    const probs = dataset.pmfByYear[c.year]
    if (!probs) continue
    for (let k = 0; k < probs.length; k++) {
      lines.push(`${c.year},${c.cohort},${k},${Math.round(probs[k] * nominalN)}`)
    }
  }
  return lines.join('\n') + '\n'
}
