# Install R packages required by the analysis notebooks and helper scripts.
#
# This list is derived from explicit library()/pkg:: usage in analysis/*.qmd
# and code/*.R files.
required <- c(
  "ipumsr",
  "dplyr",
  "tidyr",
  "readr",
  "MASS",
  "pscl",
  "knitr",
  "ggplot2",
  "MetBrewer",
  "scales",
  "purrr",
  "tibble",
  "broom",
  "viridisLite",
  "quarto",
  "rlang"
)

installed <- rownames(installed.packages())
missing <- setdiff(required, installed)

if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
  cat("Installed:", paste(missing, collapse = ", "), "\n")
} else {
  cat("All required packages are already installed.\n")
}

cat("Analysis package setup complete.\n")
