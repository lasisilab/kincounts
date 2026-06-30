#!/usr/bin/env Rscript

# Compatibility wrapper for the current IPUMS preprocessing workflow.
#
# The old version of this script expected a CSV extract and variables that are
# no longer part of the analysis. The source of truth now lives in
# analysis/preprocess.qmd, which reads the IPUMS fixed-width extract through its
# DDI/XML codebook and writes aggregate files under data/processed/.

if (!requireNamespace("quarto", quietly = TRUE)) {
  stop("Package 'quarto' is required. Install it with install.packages('quarto').", call. = FALSE)
}

if (!file.exists("analysis/preprocess.qmd")) {
  if (file.exists(file.path("..", "analysis", "preprocess.qmd"))) {
    setwd("..")
  } else {
    stop("Run this script from the project root or the code/ directory.", call. = FALSE)
  }
}

message("Rendering analysis/preprocess.qmd to rebuild or load aggregate IPUMS outputs.")
message("Set KINCOUNTS_REBUILD_IPUMS=true to force a rebuild from the local raw extract.")

quarto::quarto_render(input = "analysis/preprocess.qmd")
