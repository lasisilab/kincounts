# Install R packages required by the analysis notebooks, helper scripts,
# and optional legacy Shiny app.
#
# This list is derived from explicit library()/pkg:: usage in analysis/*.qmd,
# code/*.R, and app/*.R files.
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
  "shiny",
  "bslib",
  "DT",
  "rlang",
  "here"
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
