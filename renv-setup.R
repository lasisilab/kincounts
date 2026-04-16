# Install R packages required by analysis/*.qmd.
# This list is derived from explicit library()/pkg:: usage in analysis files.
required <- c(
  "dplyr",
  "tidyr",
  "readr",
  "MASS",
  "pscl",
  "knitr",
  "ggplot2",
  "purrr",
  "tibble",
  "broom",
  "viridisLite",
  "quarto"
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
