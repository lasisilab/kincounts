# Code Directory

This directory contains utility scripts and helper functions for the analyses:

- `utils.R`: Loads utility functions from the app directory
- `build_site.R`: Renders the Quarto website, then regenerates the simulator's empirical-parameter module
- `generate_empirical_data.R`: Generates `simulator/src/lib/empiricalData.js` from `output/fertility_parameters.csv` and `output/fertility_estimation.csv`, so the simulator's parameters never drift from the analysis
- `preprocess_ipums.R`: Earlier raw-to-processed script, now superseded by the inline processing in `analysis/preprocess.qmd`

Additional analysis-specific helper functions can be added here.
