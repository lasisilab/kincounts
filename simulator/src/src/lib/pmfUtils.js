import { zinbPMF, nbPMF } from './distributions.js'

export const BUTTERFLY_MAX_K = 12  // counts 0–11, then 12+ bin

/**
 * Compute ZINB fertility PMF for k = 0 … BUTTERFLY_MAX_K.
 * The final bin (index 12) accumulates all P(X ≥ 12).
 */
export function fertilityPMFArray(mu, theta, pi0) {
  const probs = []
  let cumulative = 0
  for (let k = 0; k < BUTTERFLY_MAX_K; k++) {
    const p = zinbPMF(k, mu, theta, pi0)
    probs.push(p)
    cumulative += p
  }
  probs.push(Math.max(0, 1 - cumulative))  // 12+ bin
  return probs  // length = BUTTERFLY_MAX_K + 1 = 13
}

/**
 * For a given fertility PMF array, derive the sibling PMF using the
 * size-biasing formula:  P(Y = k) = (k + 1) · P(X = k + 1) / E[X]
 *
 * fertProbs: array of length 13 (indices 0–11, then 12+ bin).
 * fertMean : E[X] — needed because the 12+ bin can't be decomposed further.
 *
 * The last sibling bin collects the remaining probability mass so the
 * distribution stays normalised.
 */
export function siblingPMFFromFertility(fertProbs, fertMean) {
  if (fertMean <= 0) return fertProbs.map(() => 0)

  const sibProbs = []
  let total = 0
  // k = 0 … BUTTERFLY_MAX_K - 1 (bins where P(X=k+1) is known)
  for (let k = 0; k < BUTTERFLY_MAX_K; k++) {
    const pNext = fertProbs[k + 1] ?? 0
    const p = (k + 1) * pNext / fertMean
    sibProbs.push(p)
    total += p
  }
  // 12+ bin: residual to ensure sum = 1
  sibProbs.push(Math.max(0, 1 - total))
  return sibProbs
}

/**
 * Parse a raw CSV string (from import … as 'raw' or fetch)
 * into a fertProbs array of length 13.
 * Expected format: header row + rows of "k,probability"
 */
export function parseFertilityCSV(csvText) {
  const lines = csvText.trim().split('\n').slice(1)  // skip header
  const map = {}
  for (const line of lines) {
    const [k, p] = line.split(',')
    map[Number(k)] = Number(p)
  }
  const probs = []
  let cumulative = 0
  for (let k = 0; k < BUTTERFLY_MAX_K; k++) {
    const p = map[k] ?? 0
    probs.push(p)
    cumulative += p
  }
  // Accumulate all remaining mass (k ≥ 12) into the final bin
  let tail = 0
  for (const [k, p] of Object.entries(map)) {
    if (Number(k) >= BUTTERFLY_MAX_K) tail += p
  }
  probs.push(tail > 0 ? tail : Math.max(0, 1 - cumulative))
  return probs
}

/**
 * Build the Recharts-ready data array for a butterfly panel.
 * Returns [{ count: '0', fert: -p, sib: q }, …, { count: '12+', … }]
 */
export function buildButterflyData(fertProbs, sibProbs) {
  return fertProbs.map((p, k) => ({
    count: k === BUTTERFLY_MAX_K ? '12+' : String(k),
    fert: -(p),    // negative → renders left
    sib:  sibProbs[k] ?? 0,
  }))
}
