/**
 * Return the value at a given quantile of a count array (unsorted input OK).
 * Uses a copy-sort approach; for large arrays this is the bottleneck but is
 * only called once after simulation.
 */
export function quantile(arr, q) {
  const sorted = Int32Array.from(arr).sort()
  const idx = Math.min(Math.floor(q * sorted.length), sorted.length - 1)
  return sorted[idx]
}

/**
 * Compute per-kin-type display caps from simulation results.
 * Takes allResults = { zinb, poisson, fixed } each containing typed arrays,
 * and returns { siblings: N, auntsUncles: N, ... } where each N is the
 * 99.9th percentile across all active models — so the chart window always
 * fits the actual data, no matter the parameters.
 */
export function computeDisplayCaps(allResults, kinKeys) {
  const caps = {}
  for (const key of kinKeys) {
    let maxPct = 0
    for (const modelResult of Object.values(allResults)) {
      const arr = modelResult?.[key]
      if (!arr || arr.length === 0) continue
      const pct = quantile(arr, 0.999)
      if (pct > maxPct) maxPct = pct
    }
    caps[key] = Math.max(maxPct, 1)
  }
  return caps
}

/**
 * Compute summary statistics from a raw count array.
 */
export function computeStats(arr) {
  const n = arr.length
  let sum = 0, zeros = 0
  for (let i = 0; i < n; i++) {
    sum += arr[i]
    if (arr[i] === 0) zeros++
  }
  const mean = sum / n

  let sqDev = 0
  for (let i = 0; i < n; i++) sqDev += (arr[i] - mean) ** 2

  return {
    mean,
    variance: sqDev / n,
    pZero: zeros / n,
  }
}

/**
 * Convert a raw count array to a PMF array for Recharts.
 * Returns [{ count: k, probability: p }, ...]
 */
export function arrayToPMF(arr) {
  if (!arr || arr.length === 0) return []
  const n = arr.length
  let max = 0
  for (let i = 0; i < n; i++) if (arr[i] > max) max = arr[i]

  const counts = new Float64Array(max + 1)
  for (let i = 0; i < n; i++) counts[arr[i]]++

  const pmf = []
  for (let k = 0; k <= max; k++) {
    pmf.push({ count: k, probability: counts[k] / n })
  }
  return pmf
}

/**
 * Merge multiple PMF arrays (keyed by kin type) into a single
 * Recharts data array with a shared count axis.
 *
 * mergePMFs({ siblings: PMF[], auntsUncles: PMF[] }, caps)
 *   → [{ count: 0, siblings: p, auntsUncles: q }, ...]
 *
 * caps: optional object mapping kin key → max bin to display.
 * Any probability mass beyond each key's cap is accumulated
 * into a final "X+" bin so the chart remains normalised.
 */
export function mergePMFs(pmfMap, caps = {}) {
  // Determine per-key caps and the global display cap
  let displayMax = 0
  const keyCaps = {}
  for (const key of Object.keys(pmfMap)) {
    const pmf = pmfMap[key]
    const datMax = pmf.length > 0 ? pmf[pmf.length - 1].count : 0
    const cap = caps[key] ?? datMax
    keyCaps[key] = cap
    if (cap > displayMax) displayMax = cap
  }

  const merged = []
  for (let k = 0; k <= displayMax; k++) {
    const entry = { count: k }
    for (const [key, pmf] of Object.entries(pmfMap)) {
      const keyCap = keyCaps[key]
      if (k < keyCap) {
        entry[key] = pmf[k] ? pmf[k].probability : 0
      } else if (k === keyCap) {
        // Accumulate all tail probability into the cap bin
        let tail = pmf[k] ? pmf[k].probability : 0
        for (let j = k + 1; j < pmf.length; j++) {
          tail += pmf[j] ? pmf[j].probability : 0
        }
        entry[key] = tail
      } else {
        entry[key] = 0
      }
    }
    merged.push(entry)
  }

  // Trim trailing bins where all series are effectively zero
  let last = merged.length - 1
  while (last > 0) {
    const row = merged[last]
    const allTiny = Object.entries(row)
      .filter(([k]) => k !== 'count')
      .every(([, v]) => v < 0.0005)
    if (allTiny) last--
    else break
  }
  return merged.slice(0, last + 1)
}
