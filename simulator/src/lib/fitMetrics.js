// Goodness-of-fit metrics for comparing fitted count distributions to data.
//
// Used so the Fertility Fit view reports fit from the data, never assumed:
//  - per-moment relative error (mean, variance) — the intuitive discriminator;
//  - Pearson χ² and AIC on the observed counts — the conventional whole-
//    distribution goodness-of-fit / model-selection statistics for comparing
//    Poisson vs negative-binomial / ZINB fits to count data.
//
// Note on scope: these whole-distribution statistics are only meaningful for the
// FERTILITY distribution, which is an observed sample. The sibling distribution
// is a deterministic size-biased transform of it (not an independent sample), so
// it is judged by moment errors only.

// Relative error of a model moment vs the empirical moment: |model − emp| / emp.
export function relativeError(modelValue, empValue) {
  if (empValue == null || empValue === 0) return null
  return Math.abs(modelValue - empValue) / empValue
}

// Pearson χ² over binned PMFs, using observed/expected COUNTS (O = N·emp,
// E = N·model). Returns the χ² statistic, or null if inputs are missing.
export function pearsonChiSq(empProbs, modelProbs, N) {
  if (!empProbs || !modelProbs || !N) return null
  let chi2 = 0
  for (let k = 0; k < empProbs.length; k++) {
    const E = N * (modelProbs[k] ?? 0)
    const O = N * (empProbs[k] ?? 0)
    if (E > 0) chi2 += (O - E) ** 2 / E
  }
  return chi2
}

// Multinomial log-likelihood of the binned data under a model PMF (up to the
// data-only multinomial constant, which cancels in AIC differences).
export function logLikelihood(empProbs, modelProbs, N) {
  if (!empProbs || !modelProbs || !N) return null
  let ll = 0
  for (let k = 0; k < empProbs.length; k++) {
    const O = N * (empProbs[k] ?? 0)
    const p = modelProbs[k] ?? 0
    if (O > 0 && p > 0) ll += O * Math.log(p)
  }
  return ll
}

// AIC = 2·(#params) − 2·logL. Lower is preferred; differences (ΔAIC) are what
// matter, and the dropped multinomial constant cancels in those differences.
export function aic(empProbs, modelProbs, N, nParams) {
  const ll = logLikelihood(empProbs, modelProbs, N)
  if (ll == null) return null
  return 2 * nParams - 2 * ll
}

// Parameter counts charged by AIC. ZINB nests Poisson, so on equidispersed data
// it can only match Poisson's likelihood while paying for 2 extra parameters —
// which is what makes ΔAIC prefer the simpler model when that is the truth.
export const ZINB_PARAMS = 3    // mu, theta, pi0
export const POISSON_PARAMS = 1 // lambda

/**
 * Whole-distribution comparison of the Poisson and ZINB fits to one cohort's
 * observed fertility counts: reduced Pearson χ² and ΔAIC for each, plus which
 * model AIC prefers. Shared by the Fertility Fit view and the example-fit check
 * script so the two cannot disagree.
 *
 * `nBins` is the number of PMF bins (0..11 plus the 12+ bin).
 */
export function compareFertilityFits({ empProbs, zinbProbs, poisProbs, N, nBins }) {
  if (!N) return null
  const chiZ = pearsonChiSq(empProbs, zinbProbs, N)
  const chiP = pearsonChiSq(empProbs, poisProbs, N)
  const aicZ = aic(empProbs, zinbProbs, N, ZINB_PARAMS)
  const aicP = aic(empProbs, poisProbs, N, POISSON_PARAMS)
  const minAIC = Math.min(aicZ, aicP)
  const dfZ = nBins - 1 - ZINB_PARAMS
  const dfP = nBins - 1 - POISSON_PARAMS
  return {
    N,
    zinb:    { redChi: chiZ / dfZ, dAIC: aicZ - minAIC, aic: aicZ },
    poisson: { redChi: chiP / dfP, dAIC: aicP - minAIC, aic: aicP },
    best: bestFit({ zinb: aicZ, poisson: aicP }),
  }
}

// Key of the smallest finite score (lower = better); null if nothing comparable.
export function bestFit(scores) {
  let best = null
  let bestVal = Infinity
  for (const [key, val] of Object.entries(scores)) {
    if (val == null) continue
    if (val < bestVal) { bestVal = val; best = key }
  }
  return best
}
