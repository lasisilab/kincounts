# Kincounts Simulator

An interactive browser tool accompanying the paper "Bridging fertility and sibling distributions." It shows how fertility distributions shape sibling and broader kin-count distributions, and lets you simulate the kinship network across three generations. Live: <https://kincounts.vercel.app/>.

The app has two tabs, both part of **Kincounts**:

## Fertility Fit

- Empirical fertility and sibling distributions side by side for the IPUMS cohorts (1950–1990), with ZINB and Poisson fits overlaid.
- **Computed goodness-of-fit**: each model shows a discrepancy = total-variation distance to the empirical PMF (`½·Σ|model−data|`); the closer model is marked. Nothing is hardcoded — the data decides which model wins.
- **Sibling framing**: a generation diagram makes explicit that the sibling distribution describes the *children* of the women observed at a census, not the women themselves. It is the size-biased transform P(Y=k) = (k+1)·P(X=k+1)/E[X]; the ZINB-induced sibling distribution is NB(μ·(θ+1)/θ, θ+1), in which π₀ vanishes.
- **Historical trends** plotted by census year, with the mothers' birth cohort marked under each year, and the two generations (mothers' fertility vs. children's sibship) labeled distinctly.
- **About this data** panel naming the source (IPUMS USA, CHBORN, women 50–59, 1891–1940 cohorts).
- A cohort validation table comparing model vs. empirical moments.

## Simulator

- Choose a fertility model — ZINB, Poisson, or Fixed.
- Per-generation parameters (focal / parent / grandparent), seeded from the active dataset's newest / middle / oldest cohorts.
- Simulates marginal distributions for Children, Siblings, Aunts & Uncles, Cousins, Nieces & Nephews, with mean / variance / P(0) and a reproducible seeded RNG.

## Use your own data

You can replace the IPUMS data with your own. Provide fertility as **raw counts** (women by number of children ever born, per year) via the "Use your own data" panel — upload a CSV or paste/type it, and download a pre-filled template to get the format right. The app computes the empirical distribution, fits Poisson and ZINB in-browser, and your data then drives both tabs. The fitted discrepancy is shown so you can judge fit quality.

## Export

Every chart has **CSV** (the plotted data) and **PNG** download buttons; the simulator also exports its distributions and summary statistics as CSV.

## Simulation model

**ZINB** — fertility ~ ZINB(μ, θ, π₀); sibling distribution NB(μ·(θ+1)/θ, θ+1).
**Poisson** — fertility and siblings ~ Poisson(μ).
**Fixed** — fertility and sibling counts set to their empirical means per generation.

Aunts & uncles = sum of two independent sibling draws (maternal + paternal). Cousins = random sum over aunts/uncles, each drawing from parent-generation fertility.

## Data provenance (no drift)

The IPUMS cohort parameters and PMFs are generated from the analysis outputs by `code/generate_empirical_data.R` (run from `code/build_site.R` after the Quarto render), so the simulator can never drift from the paper. Do not edit `src/lib/empiricalData.js` by hand.

## Tech stack

- [React 19](https://react.dev/) + [Vite](https://vite.dev/)
- [Recharts](https://recharts.org/) for charts
- [seedrandom](https://github.com/davidbau/seedrandom) for reproducible RNG

## Getting started

```bash
npm install
npm run dev
```

Then open http://localhost:5173. Production build: `npm run build` (outputs to `dist/`). Deployed on Vercel with the project Root Directory set to `simulator/`.
