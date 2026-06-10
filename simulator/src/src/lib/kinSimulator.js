import {
  initRNG,
  poissonSample,
  zinbSample,
  siblingFromZINB,
} from './distributions.js'

// Stochastic rounding: fills an Int32Array so values average to `mean`
function fillFixed(arr, mean) {
  const lo = Math.floor(mean)
  const hi = lo + 1
  const fracHi = mean - lo
  for (let i = 0; i < arr.length; i++) {
    arr[i] = (i / arr.length) < fracHi ? hi : lo
  }
}

/**
 * Simulate kin counts for nSim focal individuals.
 *
 * Models
 * ------
 * poisson – Poisson(mu) fertility; sibling dist is also Poisson(mu)
 * zinb    – ZINB(mu, theta, pi0) fertility; sibling dist is NB(mu*(θ+1)/θ, θ+1)
 * fixed   – fertility and sibling sizes treated as fixed at their respective means
 *           (fertMean and sibMean are independent inputs per generation)
 *
 * Generations
 * -----------
 * focal       – X₀: children, nieces/nephews
 * parent      – X₁: cousins; Y₁ (size-biased from X₁): siblings
 * grandparent – X₂: Y₂ (size-biased from X₂): aunts/uncles
 *
 * Each generation param object:
 *   ZINB:    { mu, theta, pi0 }
 *   Poisson: { mu }
 *   Fixed:   { fertMean, sibMean }
 */
export function simulateKin({ model, focal, parent, grandparent, nSim, seed }) {
  const rng = initRNG(seed)

  const children      = new Int32Array(nSim)
  const siblings      = new Int32Array(nSim)
  const auntsUncles   = new Int32Array(nSim)
  const cousins       = new Int32Array(nSim)
  const niecesNephews = new Int32Array(nSim)

  // ── Fixed model ────────────────────────────────────────────────────────────
  if (model === 'fixed') {
    // Fertility and sibling means are treated as independent fixed quantities.
    // fertMean is used for children / nieces-nephews / cousins.
    // sibMean  is used directly as the sibling / aunt-uncle count.
    const fixedChildren = focal.fertMean
    const fixedSibs     = parent.sibMean
    const fixedAU       = 2 * grandparent.sibMean
    const fixedNN       = fixedSibs  * focal.fertMean
    const fixedCos      = fixedAU    * parent.fertMean

    fillFixed(children,      fixedChildren)
    fillFixed(siblings,      fixedSibs)
    fillFixed(auntsUncles,   fixedAU)
    fillFixed(cousins,       fixedCos)
    fillFixed(niecesNephews, fixedNN)
    return { children, siblings, auntsUncles, cousins, niecesNephews }
  }

  // ── Helpers ─────────────────────────────────────────────────────────────────
  function sampleFert(gen) {
    if (model === 'poisson') return poissonSample(gen.mu, rng)
    return zinbSample(gen.mu, gen.theta, gen.pi0, rng)
  }

  function sampleSib(gen) {
    // ZINB sibling dist: NB(mu*(θ+1)/θ, θ+1) — π₀ disappears
    if (model === 'poisson') return poissonSample(gen.mu, rng)
    return siblingFromZINB(gen.mu, gen.theta, rng)
  }

  // ── Poisson / ZINB ──────────────────────────────────────────────────────────
  for (let i = 0; i < nSim; i++) {
    // Children: focal generation fertility
    const child = sampleFert(focal)

    // Siblings: size-biased from parent generation fertility
    const sibs = sampleSib(parent)

    // Aunts & uncles: 2 parents × size-biased grandparent fertility
    const auM = sampleSib(grandparent)
    const auP = sampleSib(grandparent)
    const au  = auM + auP

    // Nieces & nephews: each sibling draws from focal generation fertility
    let nn = 0
    for (let j = 0; j < sibs; j++) nn += sampleFert(focal)

    // Cousins: each aunt/uncle draws from parent generation fertility
    let cos = 0
    for (let j = 0; j < au; j++) cos += sampleFert(parent)

    children[i]      = child
    siblings[i]      = sibs
    auntsUncles[i]   = au
    cousins[i]       = cos
    niecesNephews[i] = nn
  }

  return { children, siblings, auntsUncles, cousins, niecesNephews }
}
