import seedrandom from 'seedrandom'

export function initRNG(seed) {
  return seedrandom(seed)
}

// --- Internal helpers ---

function normalSample(rng) {
  let u1
  do { u1 = rng() } while (u1 === 0)
  const u2 = rng()
  return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2)
}

// Marsaglia-Tsang gamma sampler (shape, scale)
function gammaSample(shape, scale, rng) {
  if (shape < 1) {
    return gammaSample(1 + shape, scale, rng) * Math.pow(rng(), 1 / shape)
  }
  const d = shape - 1 / 3
  const c = 1 / Math.sqrt(9 * d)
  while (true) {
    let x, v
    do {
      x = normalSample(rng)
      v = 1 + c * x
    } while (v <= 0)
    v = v * v * v
    const u = rng()
    if (u < 1 - 0.0331 * x * x * x * x) return d * v * scale
    if (Math.log(u) < 0.5 * x * x + d * (1 - v + Math.log(v))) return d * v * scale
  }
}

function logFactorial(n) {
  let sum = 0
  for (let i = 2; i <= n; i++) sum += Math.log(i)
  return sum
}

// Lanczos log-gamma (works for real z > 0)
function logGamma(z) {
  if (z < 0.5) {
    return Math.log(Math.PI / Math.sin(Math.PI * z)) - logGamma(1 - z)
  }
  z -= 1
  const c = [
    0.99999999999980993, 676.5203681218851, -1259.1392167224028,
    771.32342877765313, -176.61502916214059, 12.507343278686905,
    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7,
  ]
  let x = c[0]
  for (let i = 1; i < 9; i++) x += c[i] / (z + i)
  const t = z + 7.5
  return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x)
}

// --- Public API ---

// Poisson sample via Knuth inversion; normal approximation for lambda > 100
export function poissonSample(lambda, rng) {
  if (lambda === 0) return 0
  if (lambda > 100) {
    return Math.max(0, Math.round(lambda + Math.sqrt(lambda) * normalSample(rng)))
  }
  const L = Math.exp(-lambda)
  let k = 0, p = 1
  do { k++; p *= rng() } while (p > L)
  return k - 1
}

export function poissonPMF(k, lambda) {
  if (k < 0 || lambda <= 0) return 0
  return Math.exp(-lambda + k * Math.log(lambda) - logFactorial(k))
}

// Negative binomial sample: NB(mu, theta) — mean=mu, dispersion=theta
// Uses the Gamma-Poisson mixture: G~Gamma(theta, mu/theta), X|G~Poisson(G)
export function nbSample(mu, theta, rng) {
  if (mu === 0) return 0
  const G = gammaSample(theta, mu / theta, rng)
  return poissonSample(G, rng)
}

export function nbPMF(k, mu, theta) {
  if (k < 0) return 0
  const logC = logGamma(k + theta) - logFactorial(k) - logGamma(theta)
  const p = theta / (theta + mu)
  return Math.exp(logC + theta * Math.log(p) + k * Math.log(1 - p))
}

// ZINB sample: with prob pi0 returns 0; otherwise NB(mu, theta)
export function zinbSample(mu, theta, pi0, rng) {
  if (rng() < pi0) return 0
  return nbSample(mu, theta, rng)
}

export function zinbPMF(k, mu, theta, pi0) {
  if (k === 0) return pi0 + (1 - pi0) * nbPMF(0, mu, theta)
  return (1 - pi0) * nbPMF(k, mu, theta)
}

// Sibling distribution under ZINB fertility:
// NB(mu*(theta+1)/theta, theta+1)  — analytic result from the R codebase
export function siblingFromZINB(mu, theta, rng) {
  return nbSample(mu * (theta + 1) / theta, theta + 1, rng)
}
