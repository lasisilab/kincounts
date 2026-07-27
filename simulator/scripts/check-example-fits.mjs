#!/usr/bin/env node
// Known-truth check on the app's model selection.
//
// The two synthetic examples in src/lib/exampleData.js are generated from named
// distributions, so the correct answer is known in advance:
//
//   Poisson example       → AIC must select POISSON. ZINB nests Poisson, so it can
//                           only match the likelihood while paying for 2 extra
//                           parameters. If ZINB ever wins here, either the ZINB
//                           fitter is finding a spurious optimum or the AIC
//                           parameter counts are wrong.
//   Overdispersed example → AIC must select ZINB. Poisson cannot represent
//                           variance above the mean.
//
// This scores the same binned PMFs the Fertility Fit view plots (poissonPMFArray /
// fertilityPMFArray) with the same comparison function (compareFertilityFits), so
// a pass here means the UI shows the same verdict.
//
// Run: npm run check

import { buildUserDataset } from '../src/lib/importData.js'
import {
  poissonExampleCSV, overdispersedExampleCSV,
  POISSON_EXAMPLE, OVERDISPERSED_EXAMPLE,
} from '../src/lib/exampleData.js'
import { poissonPMFArray, fertilityPMFArray, BUTTERFLY_MAX_K } from '../src/lib/pmfUtils.js'
import { compareFertilityFits } from '../src/lib/fitMetrics.js'

const N_BINS = BUTTERFLY_MAX_K + 1

const CASES = [
  { spec: POISSON_EXAMPLE,       csv: poissonExampleCSV() },
  { spec: OVERDISPERSED_EXAMPLE, csv: overdispersedExampleCSV() },
]

let failures = 0

for (const { spec, csv } of CASES) {
  const built = buildUserDataset(csv, spec.name)
  if (built.error) {
    console.error(`FAIL  ${spec.name}: importer rejected its own example — ${built.error}`)
    failures++
    continue
  }

  console.log(`\n${spec.name}  (expected best fit: ${spec.expectedBest})`);
  console.log('  year   N       mean    var     var/mean  theta     pi0     chi2/df Z  chi2/df P  dAIC Z   dAIC P   selected');

  for (const cohort of built.dataset.cohorts) {
    const emp = built.dataset.pmfByYear[cohort.year]
    const gof = compareFertilityFits({
      empProbs:  emp,
      zinbProbs: fertilityPMFArray(cohort.mu, cohort.theta, cohort.pi0),
      poisProbs: poissonPMFArray(cohort.empMean),
      N:         cohort.sampleSize,
      nBins:     N_BINS,
    })

    const f = (v, d = 3) => Number(v).toFixed(d).padStart(8)
    console.log(
      `  ${cohort.year}  ${String(cohort.sampleSize).padStart(6)}` +
      `${f(cohort.empMean)}${f(cohort.empVariance)}${f(cohort.empVariance / cohort.empMean, 2)}` +
      `${f(cohort.theta, 1)}${f(cohort.pi0, 4)}` +
      `${f(gof.zinb.redChi, 2)}   ${f(gof.poisson.redChi, 2)}   ` +
      `${f(gof.zinb.dAIC, 2)} ${f(gof.poisson.dAIC, 2)}   ${gof.best}`
    )

    if (gof.best !== spec.expectedBest) {
      console.error(
        `  FAIL  ${spec.name} year ${cohort.year}: expected ${spec.expectedBest}, got ${gof.best} ` +
        `(AIC zinb=${gof.zinb.aic.toFixed(2)} poisson=${gof.poisson.aic.toFixed(2)})`
      )
      failures++
    }
  }
}

console.log('')
if (failures > 0) {
  console.error(`${failures} check(s) failed.`)
  process.exit(1)
}
console.log('All example-fit checks passed.')
