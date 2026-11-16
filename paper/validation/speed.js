//////////////////////////////////////////////////////////////////////////////////////
// Speed of one complete calculation (Supporting Information S5)                    //
//                                                                                  //
// One complete calculation = the simulations, the thresholds and the power of both //
// tests. It is timed for the five example scenarios of the application and for the //
// reference case, with the default numbers of simulations (1,000 selection, 2,000  //
// neutral) and one locus, and for the reference case with 1,000 candidate loci     //
// (200,000 neutral simulations).                                                   //
//                                                                                  //
// Each calculation is timed 5 times (seeds 1 to 5) after one warm-up run, on one   //
// processor core (Node.js runs on a single core). Run it with no other             //
// calculation running on the computer.                                             //
//                                                                                  //
// Run from the repository root:                                                    //
//   node paper/validation/speed.js > paper/validation/speed.tsv                    //
//////////////////////////////////////////////////////////////////////////////////////


///////////////////////////
// Load the dependencies //
///////////////////////////
const os = require('os');                    // To read the name of the processor
const S = require('../../sim.js');           // Simulation code of SPOT
const {REF} = require('./base.js');          // Settings of the reference case


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


/////////////////////////////////////////////////////////////////////////
// The scenarios timed: [name, settings]                               //
// Each example scenario = the settings above with some values changed //
/////////////////////////////////////////////////////////////////////////
// the example scenarios of the page (same values as in monte_carlo_precision.js)
const SC = [
  ['Bumblebees', D],
  ['Long-lived mammal', Object.assign({}, D, {system: 'diploid', genTime: 20, Nf: 2500, Nm: 2500, neByBlock: true, neBlocks: [5000, 200], p0: 0.1, s: 0.2, onset: 1880, years: [1880, 1950, 2020], n: [20, 20, 40], depthEarly: 4, depthLast: 15})],
  ['Human ancient DNA', Object.assign({}, D, {system: 'diploid', genTime: 28, Nf: 5000, Nm: 5000, p0: 0.1, s: 0.005, onset: -7000, years: [-7000, -5000, -3000, -1000, 1000, 2000], n: [50, 50, 50, 50, 50, 100], depthEarly: 2, depthLast: 20, err: 0.01, minDP: 1})],
  ['Bacterial pathogen', Object.assign({}, D, {system: 'haploid', genTime: 0.01, Nf: 100000, Nm: 100000, p0: 0.05, s: 0.0002, onset: 1850, years: [1850, 1900, 1950, 2000, 2020], n: [10, 10, 15, 20, 40], depthEarly: 10, depthLast: 30})],
  ['Large diploid population', Object.assign({}, D, {system: 'diploid', genTime: 1, Nf: 1e6, Nm: 1e6, p0: 0.05, s: 0.02, onset: 1900, years: [1900, 1950, 2000, 2020], n: [30, 30, 30, 60]})],
  ['Reference case', REF],
  ['Reference case, 1,000 candidate loci', Object.assign({}, REF, {m: 1000})],
];


///////////////////
// Computer used //
///////////////////
const cpu = os.cpus()[0].model.trim(),       // Name of the processor
  node = process.version;                    // Version of Node.js


///////////////////////////////////////////////////////////////////
// Time each scenario 5 times and report the median, minimum and //
// maximum (in seconds)                                          //
///////////////////////////////////////////////////////////////////
console.log(['scenario', 'loci', 'selection_simulations', 'neutral_simulations', 'median_s', 'min_s', 'max_s', 'cpu', 'node'].join('\t'));   // Header
for (const [lab, P] of SC) {                 // Loop over the scenarios
  S.runPower(Object.assign({}, P, {seed: 99}));            // warm-up
  const t = [];                              // The five run times
  for (let k = 1; k <= 5; k++) {
    const a = process.hrtime.bigint();       // Time before the calculation (nanoseconds)
    S.runPower(Object.assign({}, P, {seed: k}));   // One complete calculation
    t.push(Number(process.hrtime.bigint() - a) / 1e9);   // Time used, in seconds
  }
  t.sort((x, y) => x - y);                   // Sort the five times: t[2] is the median
  console.log([lab, P.m, P.reps, S.neutralReps(P), t[2].toFixed(3), t[0].toFixed(3), t[4].toFixed(3), cpu, node].join('\t'));   // One results line
}
