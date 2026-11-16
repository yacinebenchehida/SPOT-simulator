/////////////////////////////////////////////////////////////////////////////////
// Comparison with SLiM, part 4: power given by the threshold of SPOT on the   //
// populations simulated by SLiM (Supporting Information S4, Table S6)         //
//                                                                             //
// For each scenario with selection, one threshold is set from 20,000 neutral  //
// simulations of SPOT (95% quantile, as SPOT does). It is then applied to the //
// selection simulations of both programs, in each of the 20 independent runs  //
// (2,000 simulations per run). If SPOT simulates selection and drift as SLiM  //
// does, the threshold of SPOT gives the same power on both.                   //
//                                                                             //
// Run from the repository root, after slim_runs.sh and slim_runs_v2.sh:       //
//   node paper/validation/slim_power_spot_threshold.js                        //
// It writes paper/validation/slim_power_spot_threshold.tsv.                   //
/////////////////////////////////////////////////////////////////////////////////


///////////////////////////
// Load the dependencies //
///////////////////////////
const S = require('../../sim.js');           // Simulation code of SPOT
const fs = require('fs');                    // To read the SLiM outputs and write the results


////////////////////////////////////
// Design shared by all scenarios //
////////////////////////////////////
const GENS = [0, 20, 40, 90],                // Generations sampled
  NS = [20, 20, 20, 40],                     // Individuals sampled at each date
  R = 2000,                                  // Simulations per run
  RUNS = 20,                                 // Independent runs per scenario
  ZS = 1.6448536269514722;                   // Sampling-only threshold (95% quantile of the standard normal distribution)


/////////////////////////////////////////////////////////////////////////
// The six scenarios with selection                                    //
// [name of the SLiM run with selection, name of the neutral SLiM run, //
//  system, breeding females, breeding males, p0, s, dominance]        //
/////////////////////////////////////////////////////////////////////////
const SC = [
  ['diploid_N1000_s0.02', 'diploid_N1000_neutral', 'diploid', 500, 500, 0.1, 0.02, 0.5],
  ['diploid_N1000_s0.05_rec', 'diploid_N1000_neutral', 'diploid', 500, 500, 0.1, 0.05, 0],
  ['diploid_N200_s0.05', 'diploid_N200_neutral', 'diploid', 100, 100, 0.1, 0.05, 0.5],
  ['haploid_N1000_s0.02', 'haploid_N1000_neutral', 'haploid', 1000, 1000, 0.05, 0.02, 0.5],
  ['haplodiploid_N1000_s0.02', 'haplodiploid_N1000_neutral', 'haplodiploid', 500, 500, 0.1, 0.02, 0.5],
  ['haplodiploid_N2000_s0.05', 'haplodiploid_N2000_neutral', 'haplodiploid', 1000, 1000, 0.01, 0.05, 0.5],
];


///////////////////////////////////////////////////////////////////
// Function: read one SLiM output file                           //
// Columns: replicate, generation, pfs (population frequency), x //
// (copies of the allele sampled), nall (gene copies sampled)    //
// Returns one {pfs, x, n} per replicate                         //
///////////////////////////////////////////////////////////////////
// ---- reading SLiM output: columns replicate, generation, pfs, x, nall ----
function readSlim(name) {
  const reps = new Map();                    // One entry per replicate
  for (const line of fs.readFileSync(`paper/validation/slim_out/${name}.tsv`, 'utf8').trim().split('\n')) {   // Loop over the lines of the file
    const [id, g, pfs, x, n] = line.split('\t').map(Number),   // The five columns, as numbers
      k = GENS.indexOf(g);                   // Index of the sampling date
    if (!reps.has(id)) reps.set(id, {pfs: [], x: [], n: []});   // First line of this replicate: create its entry
    const r = reps.get(id);
    r.pfs[k] = pfs;
    r.x[k] = x;
    r.n[k] = n;
  }
  return [...reps.values()];
}


///////////////////////////////////////////////////////////////////
// Function: the same scenario simulated by SPOT (R simulations) //
// Returns one {pfs, x, n} per simulation, as readSlim           //
///////////////////////////////////////////////////////////////////
// ---- SPOT: same scenario with its own simulation code ----
function runSpot(system, Nf, Nm, p0, s, h, seed) {
  const P = S.prep({system, neByBlock: false, neBlocks: null, genTime: 1, Nf, Nm, noDrift: false, p0, s, h, onset: 0,
    years: GENS, n: NS, perfect: true, depthEarly: 10, depthLast: 15, err: 0.005, minDP: 3, alpha: 0.05, m: 1, sMale: 1,
    decline: 1, mig: 0, pMig: -1, kSib: 1, maleFrac: 0, reps: R, seed});   // Settings of the scenario (perfect genotypes, as in SLiM)
  const rng = S.makeRng(seed),               // Random number generator
    out = [];
  for (let i = 0; i < R; i++) {
    const t = S.simulateRep(P, s, rng, true);                     // One simulation
    out.push({pfs: GENS.map(g => t.traj[g]), x: t.x, n: t.n});   // Population frequency at the sampled generations, and the samples
  }
  return out;
}


////////////////////////////////////////////////////////
// Functions: statistics (the same for both programs) //
////////////////////////////////////////////////////////
// ---- statistics (identical for both simulators) ----
// Score Z of the trend test (x: copies of the allele per date; n: gene copies sampled per date)
function trendZ(x, n) {    // Cochran-Armitage trend statistic, as trendZ() in sim.js
  let N = 0,                                 // Total number of gene copies sampled
    X = 0,                                   // Total number of copies of the allele
    T = 0;                                   // Sum of the generation numbers, weighted by the gene copies sampled
  for (let i = 0; i < x.length; i++) {       // Loop over the sampling dates
    N += n[i];
    X += x[i];
    T += n[i] * GENS[i];
  }
  const pb = X / N,                          // Overall allele frequency in the samples
    tb = T / N;                              // Mean generation number of the sampled gene copies
  if (pb <= 0 || pb >= 1) return 0;          // Allele absent or fixed in all samples: Z = 0
  let num = 0,                               // Numerator: covariance between allele copies and time
    den = 0;                                 // Sum of squared deviations of time
  for (let i = 0; i < x.length; i++) {
    num += x[i] * (GENS[i] - tb);
    den += n[i] * (GENS[i] - tb) ** 2;
  }
  return num / Math.sqrt(pb * (1 - pb) * den);
}

// 95% quantile with SPOT's rule for ties at the threshold (empiricalThreshold in sim.js): {crit, gamma, ties}
const quant95 = z => S.empiricalThreshold(z, 0.05);

const mean = a => a.reduce((p, q) => p + q, 0) / a.length;   // Mean of a list of numbers

const sd = a => {                            // Standard deviation of a list of numbers
  const m = mean(a);
  return Math.sqrt(a.reduce((p, q) => p + (q - m) ** 2, 0) / (a.length - 1));
};

// Kolmogorov-Smirnov test for two samples (kept for reference; not used below)
function ks(a, b) {        // two-sample Kolmogorov-Smirnov statistic D and asymptotic p-value
  const A = [...a].sort((p, q) => p - q),    // First sample, sorted
    B = [...b].sort((p, q) => p - q);        // Second sample, sorted
  let i = 0,
    j = 0,
    D = 0;                                   // Largest distance between the two cumulative distributions
  while (i < A.length && j < B.length) {
    const v = Math.min(A[i], B[j]);
    while (i < A.length && A[i] <= v) i++;
    while (j < B.length && B[j] <= v) j++;
    D = Math.max(D, Math.abs(i / A.length - j / B.length));
  }
  const ne = A.length * B.length / (A.length + B.length),   // Effective sample size
    lam = (Math.sqrt(ne) + 0.12 + 0.11 / Math.sqrt(ne)) * D;
  let p = 0;                                 // Asymptotic p-value
  for (let k = 1; k <= 100; k++) p += 2 * (-1) ** (k - 1) * Math.exp(-2 * k * k * lam * lam);
  return [D, Math.min(1, Math.max(0, p))];
}


/////////////////
// Small tools //
/////////////////
const Z = rs => rs.map(r => trendZ(r.x, r.n));   // Score Z of each simulation of a list
// share of Z above a threshold: SPOT threshold {crit, gamma} (ties count with weight gamma), or a fixed number (sampling-only)
const above = (z, c) => typeof c === 'number' ? mean(z.map(v => v > c ? 1 : 0)) : mean(z.map(v => S.exceed(v, c.crit, c.gamma)));
const slimFile = (name, k) => k === 1 ? name : `power/${name}_r${k}`;   // File name of the SLiM run number k of a scenario


////////////////////////////////////////////////////////////////
// For each scenario: one threshold from SPOT, applied to the //
// selection simulations of SPOT and of SLiM                  //
////////////////////////////////////////////////////////////////
const rows = [['scenario', 'threshold_spot_20000', 'runs', 'mean_power_on_spot', 'sd_power_on_spot', 'mean_power_on_slim', 'sd_power_on_slim', 'difference', 'se_difference', 'z']];   // Header
SC.forEach(([sel, neu, system, Nf, Nm, p0, s, h], j) => {
  // Threshold: 95% quantile of 20,000 neutral simulations of SPOT (10 x 2,000)
  const big = [];
  for (let i = 0; i < 10; i++) big.push(...Z(runSpot(system, Nf, Nm, p0, 0, h, 50000 + 100 * j + i)));
  const c = quant95(big),                    // The threshold of SPOT
    ps = [],                                 // Power on the SPOT simulations, one value per run
    pl = [];                                 // Power on the SLiM simulations, one value per run

  // 20 runs: power of this threshold on both programs
  for (let k = 1; k <= RUNS; k++) {
    ps.push(above(Z(runSpot(system, Nf, Nm, p0, s, h, 10000 + 100 * j + k)), c));   // SPOT simulations with selection
    const zl = Z(readSlim(slimFile(sel, k)));                                       // SLiM simulations with selection
    if (zl.length !== R) throw new Error(sel + ' run ' + k);                        // Stop if a SLiM file is incomplete
    pl.push(above(zl, c));
  }

  // Difference between the two mean powers, and its standard error
  const d = mean(ps) - mean(pl),
    se = Math.sqrt(sd(ps) ** 2 / RUNS + sd(pl) ** 2 / RUNS);
  rows.push([sel, c.crit.toFixed(3), RUNS, mean(ps).toFixed(4), sd(ps).toFixed(4), mean(pl).toFixed(4), sd(pl).toFixed(4), d.toFixed(4), se.toFixed(4), (d / se).toFixed(2)]);
  console.error(sel, 'SPOT threshold', c.crit.toFixed(3), 'power on SPOT', mean(ps).toFixed(3), 'on SLiM', mean(pl).toFixed(3), 'z', (d / se).toFixed(2));   // Progress message
});


/////////////////////////////
// Write the results table //
/////////////////////////////
require('fs').writeFileSync('paper/validation/slim_power_spot_threshold.tsv', rows.map(r => r.join('\t')).join('\n') + '\n');
