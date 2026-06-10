// Per-kin-type display caps derived from the 1950 cohort (highest overdispersion):
// mu=2.818, var=6.689, theta≈2.05. Caps are ~99.5–99.9th percentile of the
// respective distributions under those empirical params. Counts beyond these
// thresholds are essentially unobservable in real census data.
export const KIN_DISPLAY_CAPS = {
  children:      12,   // fertility distribution of focal individual
  siblings:      18,   // NB(4.19, 3.05), 99.9th pct ≈ 18
  auntsUncles:   28,   // sum of 2 NB(4.19, 3.05), 99.9th pct ≈ 28
  cousins:       50,   // compound(AU, ZINB fertility), 99.5th pct ≈ 50
  niecesNephews: 45,   // compound(siblings, ZINB fertility), 99.5th pct ≈ 45
}

// IPUMS Census data — all values taken directly from the paper's GitHub repo:
//   lasisilab/kincounts  (output/fertility_estimation.csv)
//
// empMean, empVariance           — empirical fertility distribution
// zinbVariance                   — ZINB-fitted fertility variance
// empSiblingMean, empSiblingVariance  — empirical sibling distribution, constructed
//   by size-biased sampling: each mother with k+1 children contributes k+1 rows,
//   each with n_sib = k (filter(n_child > 0) %>% uncount(n_child) %>% mutate(n_sib = n_child - 1))
// zinbSiblingMean, zinbSiblingVariance — ZINB-induced sibling distribution NB(p, s+1)
//   stored directly from repo CSV to avoid rounding errors from analytical recomputation
//
// mu, theta, pi0: ZINB parameters from zeroinfl(n_child ~ 1 | 1, dist = "negbin")
//   mu  = NB component mean  (overall ZINB mean = (1−π₀)×μ)
//   theta = NB size parameter
//   pi0 = zero-inflation probability
//   Hardcoded values from simulate_relatives.qmd:
//     focal(1990): mu=3.213, theta=19.536, pi0=0.056
//     parent(1970): mu=2.530, theta=3.652, pi0=0.066
//     grandparent(1950): mu=2.943, theta=2.372, pi0=0.043
export const IPUMS_COHORTS = [
  {
    year: 1950, cohort: '1891–1900',
    empMean: 2.81794, empVariance: 6.68935,
    zinbVariance: 6.66793,
    empSiblingMean: 4.19175, empSiblingVariance: 8.89065,
    zinbSiblingMean: 4.18418, zinbSiblingVariance: 9.37581,
    mu: 2.943, theta: 2.372, pi0: 0.043,
  },
  {
    year: 1960, cohort: '1901–1910',
    empMean: 2.42926, empVariance: 5.51341,
    zinbVariance: 5.35936,
    empSiblingMean: 3.69882, empSiblingVariance: 8.51465,
    zinbSiblingMean: 3.63543, zinbSiblingVariance: 7.91489,
    mu: 2.458, theta: 2.088, pi0: 0.012,
  },
  {
    year: 1970, cohort: '1911–1920',
    empMean: 2.36310, empVariance: 4.61134,
    zinbVariance: 4.39412,
    empSiblingMean: 3.31448, empSiblingVariance: 7.09273,
    zinbSiblingMean: 3.22257, zinbSiblingVariance: 5.45518,
    mu: 2.530, theta: 3.652, pi0: 0.066,
  },
  {
    year: 1980, cohort: '1921–1930',
    empMean: 2.84623, empVariance: 4.76450,
    zinbVariance: 4.59852,
    empSiblingMean: 3.52018, empSiblingVariance: 6.30714,
    zinbSiblingMean: 3.46188, zinbSiblingVariance: 4.94432,
    mu: 3.034, theta: 7.084, pi0: 0.062,
  },
  {
    year: 1990, cohort: '1931–1940',
    empMean: 3.03355, empVariance: 4.15440,
    zinbVariance: 4.07644,
    empSiblingMean: 3.40302, empSiblingVariance: 5.03024,
    zinbSiblingMean: 3.37734, zinbSiblingVariance: 3.93278,
    mu: 3.213, theta: 19.536, pi0: 0.056,
  },
]
