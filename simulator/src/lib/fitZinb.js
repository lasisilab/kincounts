// Client-side fitting of a zero-inflated negative binomial (ZINB) to binned
// fertility counts. Used by the data-import feature so a user's own data drives
// the Fertility Fit and Simulator exactly as the IPUMS data does.
//
// The fit maximizes the multinomial log-likelihood over bins 0…11 and a 12+ tail
// (matching how the app bins fertility). It is initialized by method-of-moments
// and refined with Nelder–Mead in a transformed (unconstrained) parameter space:
//   mu = exp(a) > 0,  theta = exp(b) > 0,  pi0 = logistic(c) ∈ (0, 1).
// The Fertility Fit view shows the resulting total-variation discrepancy, so the
// quality of an in-browser fit is always visible rather than assumed.

import { zinbPMF } from './distributions.js'

const TAIL = 12 // bins 0..11 explicit, index 12 = "12+"

function logistic(x) { return 1 / (1 + Math.exp(-x)) }
function logit(p)    { return Math.log(p / (1 - p)) }

// Negative multinomial log-likelihood of binned counts under ZINB(mu, theta, pi0).
function negLL([a, b, c], counts) {
  const mu = Math.exp(a)
  const theta = Math.exp(b)
  const pi0 = logistic(c)
  let cum = 0
  let ll = 0
  for (let k = 0; k < TAIL; k++) {
    const p = zinbPMF(k, mu, theta, pi0)
    cum += p
    const n = counts[k] || 0
    if (n > 0) ll += n * Math.log(Math.max(p, 1e-12))
  }
  const pTail = Math.max(1 - cum, 1e-12)
  const nTail = counts[TAIL] || 0
  if (nTail > 0) ll += nTail * Math.log(pTail)
  return -ll
}

// Compact Nelder–Mead simplex minimizer for low-dimensional problems.
function nelderMead(f, x0, { maxIter = 600, tol = 1e-9, step = 0.25 } = {}) {
  const n = x0.length
  let simplex = [x0.slice()]
  for (let i = 0; i < n; i++) {
    const x = x0.slice()
    x[i] += x[i] !== 0 ? step * Math.abs(x[i]) : step
    simplex.push(x)
  }
  let fv = simplex.map(f)
  const alpha = 1, gamma = 2, rho = 0.5, sigma = 0.5

  for (let iter = 0; iter < maxIter; iter++) {
    const order = fv.map((v, i) => [v, i]).sort((p, q) => p[0] - q[0])
    simplex = order.map(o => simplex[o[1]])
    fv = order.map(o => o[0])
    if (Math.abs(fv[n] - fv[0]) < tol) break

    const cen = new Array(n).fill(0)
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) cen[j] += simplex[i][j]
    for (let j = 0; j < n; j++) cen[j] /= n

    const xr = cen.map((cc, j) => cc + alpha * (cc - simplex[n][j]))
    const fr = f(xr)

    if (fr >= fv[0] && fr < fv[n - 1]) { simplex[n] = xr; fv[n] = fr; continue }
    if (fr < fv[0]) {
      const xe = cen.map((cc, j) => cc + gamma * (xr[j] - cc))
      const fe = f(xe)
      if (fe < fr) { simplex[n] = xe; fv[n] = fe } else { simplex[n] = xr; fv[n] = fr }
      continue
    }
    const xc = cen.map((cc, j) => cc + rho * (simplex[n][j] - cc))
    const fc = f(xc)
    if (fc < fv[n]) { simplex[n] = xc; fv[n] = fc; continue }
    for (let i = 1; i <= n; i++) {
      for (let j = 0; j < n; j++) simplex[i][j] = simplex[0][j] + sigma * (simplex[i][j] - simplex[0][j])
      fv[i] = f(simplex[i])
    }
  }
  let bi = 0
  for (let i = 1; i < fv.length; i++) if (fv[i] < fv[bi]) bi = i
  return simplex[bi]
}

/**
 * Fit ZINB to a length-13 array of binned counts (index 0..11, index 12 = "12+").
 * `mean`/`variance` are the exact empirical moments (computed from raw counts,
 * so they may use child counts above 12). Returns { mu, theta, pi0 }.
 */
export function fitZINB(counts, mean, variance) {
  const N = counts.reduce((a, b) => a + (b || 0), 0)
  if (N === 0 || mean <= 0) return { mu: Math.max(mean, 0.01), theta: 1000, pi0: 0 }

  // ── Method-of-moments initialization ──
  // NB without inflation: if var > mean, theta0 = mean² / (var − mean).
  const theta0 = variance > mean ? Math.max((mean * mean) / (variance - mean), 0.1) : 50
  const mu0 = mean
  // Excess-zeros estimate for pi0: observed P(0) beyond what NB(mu0, theta0) predicts.
  const p0obs = (counts[0] || 0) / N
  const nbP0 = zinbPMF(0, mu0, theta0, 0)
  let pi00 = nbP0 < 1 ? (p0obs - nbP0) / (1 - nbP0) : 0
  pi00 = Math.min(Math.max(pi00, 1e-3), 0.6)

  const x = nelderMead(
    p => negLL(p, counts),
    [Math.log(mu0), Math.log(theta0), logit(pi00)]
  )

  return {
    mu: Math.exp(x[0]),
    theta: Math.min(Math.exp(x[1]), 1000), // cap to keep near-Poisson fits finite
    pi0: logistic(x[2]),
  }
}
