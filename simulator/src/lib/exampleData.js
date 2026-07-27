// Synthetic example datasets for the data importer.
//
// Every count in this file is GENERATED from a named distribution with stated
// parameters — nothing is hand-typed and nothing is census data. That matters
// twice over: the numbers can be shown in the UI without any risk of being read
// as an empirical finding, and because the generating distribution is known, the
// examples double as a check on the app's own model selection.
//
// The two examples are deliberately a contrasting pair:
//
//   Poisson       — equidispersed (variance = mean). ZINB can only match it, and
//                   pays 2 extra parameters for the privilege, so AIC should
//                   prefer POISSON.
//   Overdispersed — ZINB with real dispersion and zero-inflation. Poisson cannot
//                   represent variance > mean at all, so AIC should prefer ZINB.
//
// If either example ever selects the other model, the fit machinery has
// regressed. scripts/check-example-fits.mjs asserts exactly that.

import { poissonPMF, zinbPMF } from './distributions.js'
import { BUTTERFLY_MAX_K } from './pmfUtils.js'

const TAIL = BUTTERFLY_MAX_K // bins 0..11 explicit, index 12 = "12+"

// Bin a PMF into integer counts over 0..11 plus a 12+ bin, the same way
// aggregateByYear() bins observed data: the final bin carries P(X ≥ 12).
// Rounding drift is absorbed by the modal bin so the counts sum to exactly N.
function countsFromPMF(pmfAt, N) {
  const probs = []
  let cum = 0
  for (let k = 0; k < TAIL; k++) {
    const p = pmfAt(k)
    probs.push(p)
    cum += p
  }
  probs.push(Math.max(0, 1 - cum))

  const counts = probs.map(p => Math.round(N * p))
  const drift = N - counts.reduce((a, b) => a + b, 0)
  if (drift !== 0) {
    let mode = 0
    for (let k = 1; k < counts.length; k++) if (counts[k] > counts[mode]) mode = k
    counts[mode] += drift
  }
  return counts
}

// Serialize cohorts to the raw-counts CSV the importer accepts.
// Labels must not contain commas — the parser splits on them.
function toCSV(cohorts) {
  const lines = ['year,label,children,count']
  for (const { year, label, counts } of cohorts) {
    for (let k = 0; k < counts.length; k++) {
      lines.push(`${year},${label},${k},${counts[k]}`)
    }
  }
  return lines.join('\n') + '\n'
}

const EXAMPLE_N = 20000

// ── True Poisson: AIC should select Poisson ──
// Two cohorts so the multi-year CSV shape is visible; both equidispersed, so the
// expected selection holds whichever year you look at. The cohort label states
// the λ it was generated from, so the fitted mean is checkable on sight.
export const POISSON_EXAMPLE = {
  name: 'Synthetic Poisson',
  expectedBest: 'poisson',
  N: EXAMPLE_N,
  cohorts: [
    { year: 2000, lambda: 1.8 },
    { year: 2010, lambda: 2.2 },
  ],
}

export function poissonExampleCSV() {
  return toCSV(POISSON_EXAMPLE.cohorts.map(({ year, lambda }) => ({
    year,
    label: `λ = ${lambda.toFixed(1)}`,
    counts: countsFromPMF(k => poissonPMF(k, lambda), POISSON_EXAMPLE.N),
  })))
}

// ── Overdispersed ZINB: AIC should select ZINB ──
export const OVERDISPERSED_EXAMPLE = {
  name: 'Synthetic overdispersed',
  expectedBest: 'zinb',
  N: EXAMPLE_N,
  cohorts: [
    { year: 2000, mu: 3.0, theta: 2.5, pi0: 0.08 },
    { year: 2010, mu: 2.6, theta: 3.0, pi0: 0.12 },
  ],
}

export function overdispersedExampleCSV() {
  return toCSV(OVERDISPERSED_EXAMPLE.cohorts.map(({ year, mu, theta, pi0 }) => ({
    year,
    label: `θ = ${theta.toFixed(1)}`,
    counts: countsFromPMF(k => zinbPMF(k, mu, theta, pi0), OVERDISPERSED_EXAMPLE.N),
  })))
}

// ── Fill-in template ──
// A structurally valid, mildly overdispersed shape whose only job is to show the
// required columns and row layout. Cohort labels are obvious placeholders, and
// the numbers are synthetic so nobody mistakes the template for source data.
const TEMPLATE_N = 10000

export function makeTemplateCSV() {
  return toCSV([
    { year: 2000, label: 'Cohort A', counts: countsFromPMF(k => zinbPMF(k, 2.4, 6, 0.05), TEMPLATE_N) },
    { year: 2010, label: 'Cohort B', counts: countsFromPMF(k => zinbPMF(k, 2.0, 8, 0.05), TEMPLATE_N) },
  ])
}
