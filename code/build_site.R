# Build and preview the Quarto website
if (!requireNamespace("quarto", quietly = TRUE)) {
	stop("Package 'quarto' is required. Install with install.packages('quarto').")
}

# Render project defined by analysis/_quarto.yml
quarto::quarto_render("analysis")

# Regenerate the simulator's empirical-parameter module from the fresh model
# outputs so it cannot drift from the analysis.
source(file.path("code", "generate_empirical_data.R"))

# Uncomment to launch a live preview server
# quarto::quarto_preview("analysis")
