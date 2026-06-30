# Kincounts

Kincounts is a research project for modeling how fertility distributions shape sibling and broader kin-count distributions.

This repository includes:

- A Quarto website with methods, derivations, and analysis outputs (rendered to [docs](docs))
- An interactive **Kincounts simulator** (React + Vite), deployed on Vercel, with two parts: an empirical **Fertility Fit** view and a **Simulator** for kin counts
- A data folder that keeps restricted raw IPUMS extracts separate from shareable processed summaries
- Output folders with generated files
- A legacy R/Shiny app (kept for local use; the React simulator is the primary interactive tool)

## Project layout

- [_quarto.yml](_quarto.yml): Quarto website configuration
- [analysis](analysis): analysis pages and supporting files (including the IPUMS preprocessing methods page)
- [simulator](simulator): React + Vite interactive simulator (deployed to Vercel)
- [code](code): utility scripts, including `generate_empirical_data.R`, which syncs the simulator's empirical parameters from the analysis outputs
- [app](app): legacy R/Shiny app (local use only)
- [data](data): raw and processed data; `data/raw/` is local-only, while `data/processed/` contains shareable aggregate summaries
- [output](output): generated outputs from analyses/simulations
- [docs](docs): rendered site output
- [renv-setup.R](renv-setup.R): helper script to install analysis package dependencies

## Prerequisites

- R (recommended recent version)
- Quarto CLI installed and available in PATH
- Node.js 20+ (for the simulator)

## Setup

Run the analysis dependency setup script:

```bash
Rscript renv-setup.R
```

## Data

The raw IPUMS USA extract is not committed because IPUMS USA restricts redistribution of its data without permission. If you have IPUMS access, place the extract files in `data/raw/`. The processed aggregate summaries in `data/processed/` are produced by rendering the analysis site (the preprocessing now lives inline in [analysis/preprocess.qmd](analysis/preprocess.qmd); the earlier `code/preprocess_ipums.R` script is superseded):

```bash
quarto render
```

The repository should include only aggregate processed files in `data/processed/` so the analyses can be reproduced without redistributing record-level IPUMS data.

## Build and preview the website

Render the full Quarto site:

```bash
quarto render
```

Rendering also regenerates the simulator's empirical-parameter module (`simulator/src/lib/empiricalData.js`) from the fresh model outputs, via [code/build_site.R](code/build_site.R), so the simulator can never drift from the analysis.

Start live preview:

```bash
quarto preview
```

Rendered output is written to [docs](docs).

## Run the simulator

The interactive simulator lives in [simulator](simulator) and is deployed to Vercel.

Local development:

```bash
cd simulator
npm install
npm run dev
```

Production build:

```bash
cd simulator
npm run build   # outputs to simulator/dist
```

Vercel is configured with the project **Root Directory** set to `simulator` (so it picks up `simulator/vercel.json`, runs `npm install`, then `npm run build`).

Live app: <https://kincounts.vercel.app/>

## Run the legacy Shiny app (optional)

The React simulator is the primary interactive tool. An earlier R/Shiny app is kept for local use. To run it, install its dependencies and launch from R:

```r
install.packages(c("shiny", "bslib", "DT"), repos = "https://cloud.r-project.org")
shiny::runApp("app/app.R")
```

## Notes

- Site navigation and global rendering options are configured in [_quarto.yml](_quarto.yml).
- Analysis pages include both methodological derivations and simulation/model-fitting code.
- The repository may include generated artifacts in [docs](docs) and [output](output), depending on your workflow.
