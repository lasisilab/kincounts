# Kincounts

Kincounts is a research project for modeling how fertility distributions shape sibling and broader kin-count distributions.

This repository includes:

- A Quarto website with methods, derivations, and analysis outputs
- A Shiny app for interactive parameter estimation and simulation
- A data folder that keeps restricted raw IPUMS extracts separate from shareable processed summaries
- Output folders with generated files

## Project layout

- [_quarto.yml](_quarto.yml): Quarto website configuration
- [analysis](analysis): analysis pages and supporting files
- [app](app): Shiny app code
- [data](data): raw and processed data; `data/raw/` is local-only, while `data/processed/` contains shareable aggregate summaries
- [output](output): generated outputs from analyses/simulations
- [docs](docs): rendered site output
- [renv-setup.R](renv-setup.R): helper script to install analysis package dependencies

## Prerequisites

- R (recommended recent version)
- Quarto CLI installed and available in PATH

Optional for interactive app use:

- R packages used by the app: shiny, bslib, DT

## Setup

Run the analysis dependency setup script:

```bash
Rscript renv-setup.R
```

If you also plan to run the Shiny app, install app dependencies:

```r
install.packages(c("shiny", "bslib", "DT"), repos = "https://cloud.r-project.org")
```

## Data

The raw IPUMS USA extract is not committed because IPUMS USA restricts redistribution of its data without permission. If you have IPUMS access, place the extract files in `data/raw/` and regenerate the processed summaries:

```bash
Rscript code/preprocess_ipums.R
```

The repository should include only aggregate processed files in `data/processed/` so the analyses can be reproduced without redistributing record-level IPUMS data.

## Build and preview the website

Render the full Quarto site:

```bash
quarto render
```

Start live preview:

```bash
quarto preview
```

Rendered output is written to [docs](docs).

## Run the Shiny app

From R:

```r
shiny::runApp("app/app.R")
```

## Notes

- Site navigation and global rendering options are configured in [_quarto.yml](_quarto.yml).
- Analysis pages include both methodological derivations and simulation/model-fitting code.
- The repository may include generated artifacts in [docs](docs) and [output](output), depending on your workflow.
