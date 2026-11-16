//////////////////////////////////////////////////////////////////////////////////////////
// Precision of the power estimates (Supporting Information S5, Table S7, Figure S3)    //
//                                                                                      //
// Power is estimated from random simulations, so it varies a little between runs.      //
// This script measures how much: for the five example scenarios of the                 //
// application, the calculation is repeated with 50 different random seeds, with        //
// 500, 1,000 (default), 2,000 and 5,000 selection simulations (and twice as many       //
// neutral simulations).                                                                //
//                                                                                      //
// Run from the repository root:                                                        //
//   node paper/validation/monte_carlo_precision.js > paper/validation/mc_precision.tsv //
// It also writes paper/validation/mc_precision.json (all the power values).            //
//////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////
// Load the dependencies //
///////////////////////////
const S = require('../../sim.js');           // Simulation code of SPOT
const fs = require('fs');                    // To write the JSON file


//////////////////////////////////////////////////////////////////////
// Settings the application opens with for the "Bumblebees" example //
//////////////////////////////////////////////////////////////////////
const D = {system: 'haplodiploid',           // Inheritance system
  neByBlock: false, neBlocks: null,          // Ne is not set separately for each interval
  genTime: 1,                                // Years per generation
  Nf: 1000, Nm: 1000,                        // Breeding females and males per generation
  noDrift: false,                            // Drift is simulated
  p0: 0.01,                                  // Starting frequency of the favoured allele
  s: 0.05,                                   // Selection coefficient
  h: 0.5,                                    // Dominance
  onset: 1935,                               // Year selection starts
  years: [1935, 1955, 1975, 2025],           // Sampling years
  n: [20, 20, 20, 40],                       // Individuals sampled at each date
  perfect: false,                            // Genotypes are called from sequencing reads
  depthEarly: 10, depthLast: 15,             // Mean sequencing depth (earlier dates, last date)
  err: 0.005,                                // Error rate per read
  minDP: 3,                                  // Minimum number of reads to call a genotype
  alpha: 0.05,                               // Significance level
  m: 1,                                      // One candidate locus
  sMale: 1, decline: 1, mig: 0, pMig: -1,    // No sex difference in selection, constant Ne, no migration
  kSib: 1, maleFrac: 0,                      // Unrelated individuals, females only
  reps: 1000,                                // Selection simulations
  seed: 1};                                  // Random seed


////////////////////////////////////////////////////////////
// The five example scenarios: [name, settings]           //
// Each one = the settings above with some values changed //
////////////////////////////////////////////////////////////
// the example scenarios of the page (same values as PRESETS in parts/ui.js)
const SC = [
  ['Bumblebees', D],
  ['Long-lived mammal', Object.assign({}, D, {system: 'diploid', genTime: 20, Nf: 2500, Nm: 2500, neByBlock: true, neBlocks: [5000, 200], p0: 0.1, s: 0.2, onset: 1880, years: [1880, 1950, 2020], n: [20, 20, 40], depthEarly: 4, depthLast: 15})],
  ['Human ancient DNA', Object.assign({}, D, {system: 'diploid', genTime: 28, Nf: 5000, Nm: 5000, p0: 0.1, s: 0.005, onset: -7000, years: [-7000, -5000, -3000, -1000, 1000, 2000], n: [50, 50, 50, 50, 50, 100], depthEarly: 2, depthLast: 20, err: 0.01, minDP: 1})],
  ['Bacterial pathogen', Object.assign({}, D, {system: 'haploid', genTime: 0.01, Nf: 100000, Nm: 100000, p0: 0.05, s: 0.0002, onset: 1850, years: [1850, 1900, 1950, 2000, 2020], n: [10, 10, 15, 20, 40], depthEarly: 10, depthLast: 30})],
  ['Large diploid population', Object.assign({}, D, {system: 'diploid', genTime: 1, Nf: 1e6, Nm: 1e6, p0: 0.05, s: 0.02, onset: 1900, years: [1900, 1950, 2000, 2020], n: [30, 30, 30, 60]})],
];


/////////////////////////////////
// Fixed values and small tool //
/////////////////////////////////
const SEEDS = 50,                            // Number of repeated runs (different seeds)
  REPS = [500, 1000, 2000, 5000];            // Numbers of selection simulations
const sd = a => {                            // Function: standard deviation of a list of numbers
  const m = a.reduce((x, y) => x + y, 0) / a.length;                                  // Mean
  return Math.sqrt(a.reduce((x, y) => x + (y - m) ** 2, 0) / (a.length - 1));         // Square root of the variance
};
const all = {};                              // All the power values, kept for the JSON file


////////////////////////////////////////////////////////////////////
// For each scenario and number of simulations: 50 runs, then the //
// mean and the standard deviation of the 50 power estimates      //
////////////////////////////////////////////////////////////////////
console.log(['scenario', 'replicates', 'statistic', 'mean_power', 'sd_power', 'sd_binomial_only', 'sd_ratio', 'min', 'max', 'seeds'].join('\t'));   // Header
const t0 = Date.now();                       // Starting time, to print the time used
for (const [lab, P] of SC) {                 // Loop over the scenarios
  all[lab] = {};
  for (const R of REPS) {                    // Loop over the numbers of selection simulations
    const v = {trend: [], two: []};          // Power of the trend test and of the first-versus-last test, one value per run
    for (let i = 0; i < SEEDS; i++) {        // 50 runs with different seeds
      const r = S.runPower(Object.assign({}, P, {reps: R, seed: 5000 + 97 * i}));   // One SPOT run
      v.trend.push(r.trend.powerEmp);
      v.two.push(r.two.powerEmp);
    }
    all[lab][R] = v;
    for (const k of ['trend', 'two']) {      // For both tests
      const a = v[k],                        // The 50 power estimates
        m = a.reduce((x, y) => x + y, 0) / a.length,   // Their mean
        b = Math.sqrt(m * (1 - m) / R);      // Standard deviation expected if the threshold did not vary (binomial only)
      console.log([lab, R, k, m.toFixed(4), sd(a).toFixed(4), b.toFixed(4), (sd(a) / b).toFixed(2), Math.min(...a).toFixed(3), Math.max(...a).toFixed(3), SEEDS].join('\t'));   // One results line
    }
    console.error(lab, R, ((Date.now() - t0) / 1000).toFixed(0) + 's');   // Progress message with the time used
  }
}


///////////////////////////////
// Save all the power values //
///////////////////////////////
fs.writeFileSync('paper/validation/mc_precision.json', JSON.stringify(all));
