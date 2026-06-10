# Kin Counts Simulator

An interactive browser tool accompanying the paper "Bridging fertility and sibling distributions." It visualizes the relationship between fertility and sibling distributions under different models, and simulates the full kinship network across three generations.

## Sections

### Fertility Fit
- Displays empirical fertility and sibling distributions side by side for five IPUMS cohorts (1950–1990)
- Overlays ZINB and Poisson model fits on both charts
- Sibling distribution constructed via size-biasing: P(Y=k) = (k+1)·P(X=k+1)/E[X]
- ZINB-fitted sibling distribution: NB(μ·(θ+1)/θ, θ+1) — zero-inflation π₀ vanishes
- Stats panel reports mean and variance for empirical, ZINB, and Poisson across both distributions
- Cohort validation table reproduces Table 2 from the paper across all five cohorts

### Kinship Counts Simulator
- Choose a fertility model: ZINB, Poisson, or Fixed
- Set parameters independently for three generations: focal (1990), parent (1970), grandparent (1950)
- Simulates 100,000 individuals and produces full marginal distributions for:
  - Children, Siblings, Aunts & Uncles, Cousins, Nieces & Nephews
- Summary statistics: mean, variance, and P(0) per kin type
- Reproducible results via seeded RNG

## Simulation Model

**ZINB** — fertility drawn from ZINB(μ, θ, π₀); sibling distribution is NB(μ·(θ+1)/θ, θ+1)

**Poisson** — fertility and siblings drawn from Poisson(μ)

**Fixed** — fertility and sibling counts set to independent empirical means per generation

Aunts & uncles = sum of two independent sibling draws (maternal + paternal). Cousins = random sum over aunts/uncles, each drawing from parent-generation fertility.

## Default Parameters

Defaults match fitted values from the paper (IPUMS USA, women aged 50–59):

| Generation | Model | μ | θ | π₀ |
|------------|-------|---|---|-----|
| Focal (1990) | ZINB | 3.213 | 19.536 | 0.056 |
| Parent (1970) | ZINB | 2.530 | 3.652 | 0.066 |
| Grandparent (1950) | ZINB | 2.943 | 2.372 | 0.043 |

## Tech Stack

- [React 19](https://react.dev/) + [Vite](https://vite.dev/)
- [Recharts](https://recharts.org/) for charts
- [seedrandom](https://github.com/davidbau/seedrandom) for reproducible RNG

## Getting Started

```bash
npm install
npm run dev
```

Then open http://localhost:5173 in your browser.
