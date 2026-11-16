//////////////////////////////////////////////////////////////////////////////////////
// Rescaling: false-positive rate of the drift-aware threshold with rescaling (S6)  //
//                                                                                  //
// Question: when SPOT sets its threshold from rescaled neutral simulations, is     //
// the false-positive rate still the intended 5%?                                   //
//                                                                                  //
// For each design, SPOT sets its threshold from rescaled neutral simulations       //
// (20 independent runs, with the default 1,000 selection and 2,000 neutral         //
// simulations each). The realised false-positive rate of each threshold is then    //
// measured on 50,000 independent neutral simulations run WITHOUT rescaling (every  //
// generation simulated). The sampling-only threshold is evaluated on the same set. //
//                                                                                  //
// Settings: four sampling dates equally spaced over G generations, 20/20/20/40     //
// individuals, default sequencing settings, p0 = 0.1, breeding females = breeding  //
// males = 10,000.                                                                  //
//                                                                                  //
// Run from the repository root:                                                    //
//   node paper/validation/rescale_fpr.js > paper/validation/rescale_fpr.tsv        //
//////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////
// Load the simulation code //
//////////////////////////////
const S = require('../../sim.js');


////////////////////////////////////
// Settings shared by all designs //
////////////////////////////////////
const base = {genTime: 1,            // Years per generation
  Nf: 1e4, Nm: 1e4,                  // Breeding females and males per generation
  noDrift: false,                    // Drift is simulated
  p0: 0.1,                           // Starting frequency of the allele
  s: 0,                              // No selection (only neutral simulations matter here)
  h: 0.5,                            // Dominance
  onset: 0,                          // Start of selection (not used, s = 0)
  n: [20, 20, 20, 40],               // Individuals sampled at each date
  perfect: false,                    // Genotypes are called from sequencing reads
  depthEarly: 10, depthLast: 15,     // Mean sequencing depth (earlier dates, last date)
  err: 0.005,                        // Error rate per read
  minDP: 3,                          // Minimum number of reads to call a genotype
  alpha: 0.05,                       // Significance level
  m: 1,                              // One candidate locus
  sMale: 1, decline: 1, mig: 0, pMig: -1,   // No sex difference in selection, constant Ne, no migration
  kSib: 1,                           // Unrelated individuals
  maleFrac: 0,                       // Females only
  reps: 1000,                        // Selection simulations per run (default)
  seed: 1};                          // Random seed
const N_THRESHOLDS = 20,             // Number of independent SPOT runs (one threshold each)
  N_TEST = 50000;                    // Number of independent neutral simulations without rescaling


////////////////////////////////////////////////////////////////////////
// Function: score Z of the trend test                                //
// x: copies of the allele per date; n: gene copies sampled per date; //
// t: generation number of each date; idx: dates used by the test     //
////////////////////////////////////////////////////////////////////////
function trendZ(x, n, t, idx) {
  let N = 0,                         // Total number of gene copies sampled
    X = 0,                           // Total number of copies of the allele
    T = 0;                           // Sum of the generation numbers, weighted by the gene copies sampled
  for (const i of idx) {             // Loop over the sampling dates
    N += n[i];
    X += x[i];
    T += n[i] * t[i];
  }
  if (N === 0) return 0;             // No gene copy sampled: Z = 0
  const pb = X / N,                  // Overall allele frequency in the samples
    tb = T / N;                      // Mean generation number of the sampled gene copies
  if (pb <= 0 || pb >= 1) return 0;  // Allele absent or fixed in all samples: Z = 0
  let a = 0,                         // Numerator: covariance between allele copies and time
    b = 0;                           // Sum of squared deviations of time
  for (const i of idx) {             // Loop over the sampling dates
    a += x[i] * (t[i] - tb);
    b += n[i] * (t[i] - tb) ** 2;
  }
  return b <= 0 ? 0 : a / Math.sqrt(pb * (1 - pb) * b);   // Cochran-Armitage trend score
}


//////////////////////////////////////////////////////////////////////
// For each design: thresholds from rescaled simulations, tested on //
// neutral simulations without rescaling                            //
//////////////////////////////////////////////////////////////////////
console.log(['G', 'system', 'Ne', 'lambda', 'fpr_mean', 'fpr_sd', 'fpr_min', 'fpr_max', 'fpr_sampling_only', 'threshold_mean', 'n_test'].join('\t'));   // Header
const t0 = Date.now();                                                 // Starting time, to print the time used
for (const G of [5000, 20000, 60000])                                  // Loop over the numbers of generations
  for (const system of ['haploid', 'diploid', 'haplodiploid']) {       // Loop over the inheritance systems
  const P = Object.assign({}, base, {system, years: [0, Math.round(G / 3), Math.round(2 * G / 3), G]});   // Settings of this design (four equally spaced dates)
  const lambda = S.rescale(P).lambda;                                  // Scaling factor chosen by SPOT

  // 20 thresholds, each from one SPOT run with rescaling
  const crit = [];
  for (let i = 0; i < N_THRESHOLDS; i++)
    crit.push(S.runPower(Object.assign({}, P, {seed: 1000 + i})).trend);   // {crit, gamma}

  // 50,000 independent neutral simulations WITHOUT rescaling, and their score Z
  const Q = S.prep(Object.assign({}, P, {noRescale: true})),           // Settings without rescaling, in the form used by the simulation code
    R = S.makeRng(99991),                                              // Independent random number generator
    t = S.gensOf(Q),                                                   // Generation number of each sampling date
    ALL = t.map((_, i) => i);                                          // Indices of all sampling dates
  const z = new Float64Array(N_TEST);                                  // Score Z of each neutral simulation
  for (let i = 0; i < N_TEST; i++) {
    const o = S.simulateRep(Q, 0, R, false);                           // One neutral simulation
    z[i] = trendZ(o.x, o.n, t, ALL);                                   // Its trend score
  }

  // False-positive rate of each of the 20 thresholds
  const f = crit.map(c => {
    let k = 0;
    for (const v of z) k += S.exceed(v, c.crit, c.gamma);
    return k / N_TEST;
  });   // ties count with weight gamma
  const m = f.reduce((a, b) => a + b, 0) / f.length,                   // Mean of the 20 false-positive rates
    sd = Math.sqrt(f.reduce((a, b) => a + (b - m) ** 2, 0) / (f.length - 1));   // Their standard deviation

  // False-positive rate of the sampling-only threshold on the same simulations
  const zc = S.qnorm(0.95);                                            // Sampling-only threshold (1.645)
  let naive = 0;
  for (const v of z) if (v > zc) naive++;                              // Count the neutral simulations above it

  // One results line per design
  console.log([G, system, Math.round(S.neFor(P)), lambda, m.toFixed(4), sd.toFixed(4), Math.min(...f).toFixed(4), Math.max(...f).toFixed(4),
    (naive / N_TEST).toFixed(4), (crit.reduce((a, b) => a + b.crit, 0) / crit.length).toFixed(3), N_TEST].join('\t'));
  console.error(G, system, 'lambda', lambda, 'FPR', m.toFixed(4), ((Date.now() - t0) / 1000).toFixed(0) + 's');   // Progress message with the time used
}
