#!/usr/bin/env Rscript
# Generate simulator/src/lib/empiricalData.js from the fertility output CSVs, so
# the simulator's empirical parameters can never drift from the analysis.
#
# Source of truth:
#   output/fertility_parameters.csv  (ZINB mu/size/pi0, written by fert_model.qmd)
#   output/fertility_estimation.csv  (empirical & ZINB means/variances)
#
# Run after fert_model.qmd has been rendered (it writes those CSVs), e.g. via
# code/build_site.R. Do not edit empiricalData.js by hand.

suppressMessages({
  library(readr)
  library(dplyr)
})

root <- if (dir.exists("output")) "." else ".."
out_dir <- file.path(root, "output")
js_path <- file.path(root, "simulator", "src", "lib", "empiricalData.js")

params <- read_csv(file.path(out_dir, "fertility_parameters.csv"), show_col_types = FALSE)
estimation <- read_csv(file.path(out_dir, "fertility_estimation.csv"), show_col_types = FALSE)

# Per-cohort sample sizes (unweighted) so the simulator can report real
# goodness-of-fit statistics (Pearson chi-square, AIC) on the observed counts.
sizes <- read_csv(
  file.path(root, "data", "processed", "ipums_cohort_sample_sizes.csv"),
  show_col_types = FALSE
)

# Cohort display labels are nominal descriptors (they do not drift with the data),
# so they stay as a fixed lookup rather than being read from a CSV.
cohort_labels <- c(
  "1950" = "1891–1900",
  "1960" = "1901–1910",
  "1970" = "1911–1920",
  "1980" = "1921–1930",
  "1990" = "1931–1940"
)

fert <- estimation %>% filter(distribution == "fertility")
sib  <- estimation %>% filter(distribution == "sibling")
zinb <- params %>% filter(Model == "ZINB")

years <- sort(unique(zinb$Year))

entries <- vapply(years, function(y) {
  f <- fert[fert$Year == y, ]
  s <- sib[sib$Year == y, ]
  z <- zinb[zinb$Year == y, ]
  label <- cohort_labels[[as.character(y)]]
  if (is.null(label)) stop("No cohort label defined for year ", y)
  n_women <- sizes$n_women[match(y, sizes$Year)]
  if (is.na(n_women)) stop("No sample size found for year ", y)
  sprintf(
    paste(
      "  {",
      "    year: %d, cohort: '%s',",
      "    empMean: %.5f, empVariance: %.5f,",
      "    zinbVariance: %.5f,",
      "    empSiblingMean: %.5f, empSiblingVariance: %.5f,",
      "    zinbSiblingMean: %.5f, zinbSiblingVariance: %.5f,",
      "    mu: %.3f, theta: %.3f, pi0: %.3f,",
      "    sampleSize: %d,",
      "  },",
      sep = "\n"
    ),
    y, label,
    f$mean_data, f$variance_data, f$variance_zinb,
    s$mean_data, s$variance_data, s$mean_zinb, s$variance_zinb,
    z$mu, z$size, z$pi0,
    as.integer(n_women)
  )
}, character(1))

header <- c(
  "// GENERATED FILE — do not edit by hand.",
  "// Produced by code/generate_empirical_data.R from the fertility output CSVs",
  "// (output/fertility_parameters.csv, output/fertility_estimation.csv). Re-run that",
  "// script, or code/build_site.R, after fert_model.qmd to refresh these values.",
  "//",
  "// Per-kin-type display caps derived from the 1950 cohort (highest overdispersion).",
  "// These are fixed ~99.5–99.9th-percentile display thresholds, not regenerated",
  "// from the CSVs; counts beyond them are essentially unobservable in census data.",
  "export const KIN_DISPLAY_CAPS = {",
  "  children:      12,",
  "  siblings:      18,",
  "  auntsUncles:   28,",
  "  cousins:       50,",
  "  niecesNephews: 45,",
  "}",
  "",
  "// IPUMS Census cohorts. All numeric fields are generated from the output CSVs:",
  "//   empMean, empVariance              — empirical fertility mean and variance",
  "//   zinbVariance                      — ZINB-fitted fertility variance",
  "//   empSiblingMean, empSiblingVariance  — empirical (size-biased) sibling distribution",
  "//   zinbSiblingMean, zinbSiblingVariance — ZINB-induced sibling distribution",
  "//   mu, theta, pi0                    — ZINB parameters (overall mean = (1−pi0)*mu)",
  "//   sampleSize                        — unweighted n women (for chi-square / AIC)",
  "export const IPUMS_COHORTS = ["
)

js <- c(header, entries, "]", "")

con <- file(js_path, open = "w", encoding = "UTF-8")
writeLines(js, con)
close(con)

# Sync the empirical PMF curves the simulator imports (the other census-data source),
# so both come from the regenerated outputs and cannot drift.
pmf_src <- Sys.glob(file.path(root, "data", "processed", "fertility_pmf_*.csv"))
sim_data <- file.path(root, "simulator", "src", "data")
if (length(pmf_src) > 0 && dir.exists(sim_data)) {
  file.copy(pmf_src, sim_data, overwrite = TRUE)
}

message("Wrote ", js_path, " and synced ", length(pmf_src), " PMF CSVs to the simulator.")
