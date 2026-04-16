# Build and preview the Quarto website
if (!requireNamespace("quarto", quietly = TRUE)) {
	stop("Package 'quarto' is required. Install with install.packages('quarto').")
}

# Render project defined by analysis/_quarto.yml
quarto::quarto_render("analysis")

# Uncomment to launch a live preview server
# quarto::quarto_preview("analysis")
