/////////////////////////////////////////////////////////////////////////////////////////////
// Reference case of the article (Figure 1; p0 = 0.1, s = 0.02), run 20 times              //
// (abstract and Section 3.4 of the main text)                                             //
//                                                                                         //
// Figure 1 shows one run with 1,000 selection simulations. This script repeats the        //
// calculation 20 times with different seeds, with 1,000 and with 5,000 selection          //
// simulations, and reports the mean, standard deviation and range over the 20 runs of:    //
//   - the power of the trend test with the drift-aware threshold,                         //
//   - the power the sampling-only threshold would report,                                 //
//   - the false-positive rate of the sampling-only threshold.                             //
//                                                                                         //
// Run from the repository root:                                                           //
//   node paper/worked_examples/headline_runs.js > paper/worked_examples/headline_runs.tsv //
/////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////
// Load the simulation code //
//////////////////////////////
const S = require('../../sim.js');


///////////////////////////////////////////////////////////
// Settings of the reference case (as in the app, Fig 1) //
///////////////////////////////////////////////////////////
const REF = {system: 'haplodiploid',         // Inheritance system
  neByBlock: false, neBlocks: null,          // Ne is not set separately for each interval
  genTime: 1,                                // Years per generation
  Nf: 1000, Nm: 1000,                        // Breeding females and males per generation (Ne = 1,500)
  noDrift: false,                            // Drift is simulated
  p0: 0.1,                                   // Starting frequency of the favoured allele
  s: 0.02,                                   // Selection coefficient
  h: 0.5,                                    // Dominance
  onset: 1935,                               // Year selection starts (the first sampling date)
  years: [1935, 1955, 1975, 2025],           // Sampling years
  n: [20, 20, 20, 40],                       // Individuals sampled at each date
  perfect: false,                            // Genotypes are called from sequencing reads
  depthEarly: 10, depthLast: 15,             // Mean depth at the earlier dates and at the last date
  err: 0.005, errLast: 0.005,                // Per-read error rate at the earlier dates and at the last date
  minDP: 3,                                  // Minimum depth to call a genotype
  alpha: 0.05, m: 1,                         // Significance level and number of loci tested
  sMale: 1,                                  // Selection in males (x s)
  decline: 1,                                // Ne at the last date relative to the first
  mig: 0, pMig: -1,                          // No migration
  kSib: 1,                                   // Unrelated individuals
  maleFrac: 0};                              // Females only


////////////////////////////////////////////////
// Functions: mean, standard deviation, range //
////////////////////////////////////////////////
const mean = v => v.reduce((a, b) => a + b, 0) / v.length;                                 // Mean of a list
const sd = v => Math.sqrt(v.reduce((a, b) => a + (b - mean(v)) ** 2, 0) / (v.length - 1)); // Standard deviation of a list
const pc = x => (100 * x).toFixed(1);                                                      // Proportion written as a percentage, one decimal


/////////////////////////////////////////////////////
// 20 runs (seeds 1 to 20) for each number of      //
// selection simulations, then one line per number //
/////////////////////////////////////////////////////
console.log(['selection_simulations', 'runs', 'measure', 'mean_percent', 'sd_percent', 'min_percent', 'max_percent'].join('\t'));   // Header of the results file
for (const reps of [1000, 5000]) {           // Default number of selection simulations, and the number recommended for final estimates
  const da = [], so = [], fp = [];           // Drift-aware power, sampling-only power, sampling-only false-positive rate
  for (let seed = 1; seed <= 20; seed++) {   // 20 runs with different seeds
    const r = S.runPower(Object.assign({}, REF, {reps, seed}));   // One complete calculation, as in the app
    da.push(r.trend.powerEmp);               // Power of the trend test, drift-aware threshold
    so.push(r.trend.powerNaive);             // Power reported by the sampling-only threshold
    fp.push(r.trend.fprNaive);               // False-positive rate of the sampling-only threshold
  }
  for (const [name, v] of [['power_drift_aware', da], ['power_sampling_only', so], ['false_positive_rate_sampling_only', fp]])   // One line per measure
    console.log([reps, 20, name, pc(mean(v)), pc(sd(v)), pc(Math.min(...v)), pc(Math.max(...v))].join('\t'));
}
