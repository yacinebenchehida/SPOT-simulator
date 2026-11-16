///////////////////////////////////////////////////////////////////////////////////////////
// Example without drift (Supporting Information S1.11)                                  //
//                                                                                       //
// Question: without drift, do the drift-aware threshold and the sampling-only threshold //
// become identical? Answer: almost, but not exactly, when the favoured allele is rare,  //
// because the score Z can then take only a few distinct values.                         //
//                                                                                       //
// Settings: diploid, no drift, perfect genotypes, unrelated individuals, p0 = 0.01,     //
// s = 0.05, sampling design of the reference case (museum insect case study).           //
//                                                                                       //
// Run from the repository root:                                                         //
//   node paper/validation/nodrift_example.js > paper/validation/nodrift_example.tsv     //
///////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////
// Load the simulation code and the settings //
///////////////////////////////////////////////
const {S, with: W} = require('./base.js');   // S: the simulation code (sim.js); W: the reference case with some settings changed
const P = W({                                // Settings of this example (everything else as in the reference case)
  system: 'diploid',                         // Diploid inheritance
  noDrift: true,                             // Drift switched off (infinite population)
  perfect: true,                             // Perfect genotypes (no sequencing step)
  p0: 0.01,                                  // Rare favoured allele at the first date
  s: 0.05,                                   // Selection coefficient
  reps: 1000,                                // 1,000 selection simulations (and 2,000 neutral ones, the default)
  seed: 1                                    // Random seed, so that the run can be reproduced
});


/////////////////////////////////////////
// One SPOT run: thresholds and powers //
/////////////////////////////////////////
const r = S.runPower(P);                     // Run SPOT once with these settings
const crit = r.trend.crit;                   // Drift-aware threshold of the trend test
const gamma = r.trend.gamma;                 // Weight given to simulations whose Z equals the threshold (ties; see sim.js)
const ties = r.trend.ties;                   // Number of neutral simulations tied with the threshold
const zc = r.zNaive;                         // Sampling-only threshold (1.645 for alpha = 0.05 and one locus)


/////////////////////////////////////////////
// Prepare independent neutral simulations //
/////////////////////////////////////////////
const prep = S.prep(S.rescale(P).P);         // Settings in the form used by the simulation code
const R = S.makeRng(777);                    // A second random number generator, independent of the SPOT run above
const t = S.gensOf(prep);                    // Generation number of each sampling date
const all = t.map((_, i) => i);              // Indices of all sampling dates (the trend test uses all of them)


///////////////////////////////////////////////
// Function: score Z of the trend test       //
// x: copies of the favoured allele per date //
// n: gene copies sampled per date           //
///////////////////////////////////////////////
const trendZ = (x, n) => {
  let N = 0;                                 // Total number of gene copies sampled
  let X = 0;                                 // Total number of copies of the favoured allele
  let T = 0;                                 // Sum of the generation numbers, weighted by the gene copies sampled
  for (const i of all) {                     // Loop over the sampling dates
    N += n[i];
    X += x[i];
    T += n[i] * t[i];
  }
  if (!N) return 0;                          // No gene copy sampled: Z = 0
  const pb = X / N;                          // Overall frequency of the favoured allele in the samples
  const tb = T / N;                          // Mean generation number of the sampled gene copies
  if (pb <= 0 || pb >= 1) return 0;          // Allele absent or fixed in all samples: Z = 0
  let a = 0;                                 // Numerator: covariance between allele copies and time
  let b = 0;                                 // Sum of squared deviations of time
  for (const i of all) {                     // Loop over the sampling dates
    a += x[i] * (t[i] - tb);
    b += n[i] * (t[i] - tb) ** 2;
  }
  return b > 0 ? a / Math.sqrt(pb * (1 - pb) * b) : 0;   // Cochran-Armitage trend score
};


/////////////////////////////////////////////////////////////////
// Count the selection simulations above or equal to threshold //
/////////////////////////////////////////////////////////////////
let sGt = 0;                                 // Selection simulations strictly above the drift-aware threshold
let sEq = 0;                                 // Selection simulations tied with the drift-aware threshold
for (const z of r.sel.trend) {               // Loop over the Z of the 1,000 selection simulations
  if (S.isTie(z, crit)) sEq++;               // Tied with the threshold
  else if (z > crit) sGt++;                  // Strictly above the threshold
}


///////////////////////////////////////////////////////////////////
// False-positive rate of each threshold, on 100,000 independent //
// neutral simulations                                           //
///////////////////////////////////////////////////////////////////
const NT = 100000;                           // Number of independent neutral simulations
let f1 = 0;                                  // Neutral simulations counted above the drift-aware threshold (ties weighted by gamma)
let f2 = 0;                                  // Neutral simulations above the sampling-only threshold
for (let i = 0; i < NT; i++) {               // Loop over the neutral simulations
  const o = S.simulateRep(prep, 0, R, false);   // One neutral simulation (selection coefficient 0)
  const z = trendZ(o.x, o.n);                // Its trend score
  f1 += S.exceed(z, crit, gamma);            // 1 if above, gamma if tied, 0 if below
  if (z > zc) f2++;                          // Above the sampling-only threshold
}


///////////////////////////////////////
// Write the results (tab-separated) //
///////////////////////////////////////
console.log([                                // Header line
  'threshold_neutral_simulations',
  'n_ties_at_threshold',
  'gamma_ties',
  'threshold_sampling_only',
  'power_neutral_simulation_threshold',
  'power_sampling_only',
  'fpr_neutral_simulation_threshold',
  'fpr_sampling_only',
  'n_independent_neutral',
  'n_selection_above',
  'n_selection_tied',
  'power_strictly_above',
  'power_ties_all_counted'
].join('\t'));
console.log([                                // Results line
  crit.toFixed(3),                           // Drift-aware threshold
  ties,                                      // Neutral simulations tied with it
  gamma.toFixed(4),                          // Weight of the ties
  zc.toFixed(3),                             // Sampling-only threshold
  r.trend.powerEmp.toFixed(4),               // Power with the drift-aware threshold
  r.trend.powerNaive.toFixed(4),             // Power with the sampling-only threshold
  (f1 / NT).toFixed(4),                      // False-positive rate of the drift-aware threshold
  (f2 / NT).toFixed(4),                      // False-positive rate of the sampling-only threshold
  NT,                                        // Number of independent neutral simulations
  sGt,                                       // Selection simulations strictly above the threshold
  sEq,                                       // Selection simulations tied with the threshold
  (sGt / P.reps).toFixed(4),                 // Power if ties are not counted
  ((sGt + sEq) / P.reps).toFixed(4)          // Power if all ties are counted
].join('\t'));
