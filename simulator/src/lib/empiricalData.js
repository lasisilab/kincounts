// GENERATED FILE — do not edit by hand.
// Produced by code/generate_empirical_data.R from the fertility output CSVs
// (output/fertility_parameters.csv, output/fertility_estimation.csv). Re-run that
// script, or code/build_site.R, after fert_model.qmd to refresh these values.
//
// Per-kin-type display caps derived from the 1950 cohort (highest overdispersion).
// These are fixed ~99.5–99.9th-percentile display thresholds, not regenerated
// from the CSVs; counts beyond them are essentially unobservable in census data.
export const KIN_DISPLAY_CAPS = {
  children:      12,
  siblings:      18,
  auntsUncles:   28,
  cousins:       50,
  niecesNephews: 45,
}

// IPUMS Census cohorts. All numeric fields are generated from the output CSVs:
//   empMean, empVariance              — empirical fertility mean and variance
//   zinbVariance                      — ZINB-fitted fertility variance
//   empSiblingMean, empSiblingVariance  — empirical (size-biased) sibling distribution
//   zinbSiblingMean, zinbSiblingVariance — ZINB-induced sibling distribution
//   mu, theta, pi0                    — ZINB parameters (overall mean = (1−pi0)*mu)
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
    zinbVariance: 5.35937,
    empSiblingMean: 3.69882, empSiblingVariance: 8.51465,
    zinbSiblingMean: 3.63543, zinbSiblingVariance: 7.91489,
    mu: 2.458, theta: 2.088, pi0: 0.012,
  },
  {
    year: 1970, cohort: '1911–1920',
    empMean: 2.37699, empVariance: 4.67520,
    zinbVariance: 4.47304,
    empSiblingMean: 3.34383, empSiblingVariance: 7.05668,
    zinbSiblingMean: 3.25880, zinbSiblingVariance: 5.57870,
    mu: 2.547, theta: 3.578, pi0: 0.067,
  },
  {
    year: 1980, cohort: '1921–1930',
    empMean: 2.85487, empVariance: 4.78726,
    zinbVariance: 4.62013,
    empSiblingMean: 3.53174, empSiblingVariance: 6.34570,
    zinbSiblingMean: 3.47320, zinbSiblingVariance: 4.99107,
    mu: 3.036, theta: 6.947, pi0: 0.060,
  },
  {
    year: 1990, cohort: '1931–1940',
    empMean: 3.04652, empVariance: 4.16000,
    zinbVariance: 4.08293,
    empSiblingMean: 3.41201, empSiblingVariance: 5.02632,
    zinbSiblingMean: 3.38672, zinbSiblingVariance: 3.93837,
    mu: 3.224, theta: 19.792, pi0: 0.055,
  },
]

