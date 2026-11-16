/////////////////////////////////////////////////////////////////////////////////////
// Agreement with theory: the values given by SPOT (Supporting Information S3)     //
//                                                                                 //
// This script produces the SPOT side of the comparison with theory. The           //
// theoretical values are computed independently by theory_check.py, which then    //
// compares the two. Three checks:                                                 //
//   (a) power without drift and with perfect genotypes (only sampling is random); //
//   (b) loss of neutral genetic diversity by drift, generation by generation;     //
//   (c) allele-frequency trajectories under selection without drift.              //
//                                                                                 //
// Run from the repository root:                                                   //
//   node paper/validation/theory_spot.js                                          //
// It writes paper/validation/theory_spot.json.                                    //
/////////////////////////////////////////////////////////////////////////////////////


///////////////////////////
// Load the dependencies //
///////////////////////////
const S = require('../../sim.js');           // Simulation code of SPOT
const fs = require('fs');                    // To write the JSON file


/////////////////////////////////////////
// Settings shared by the three checks //
/////////////////////////////////////////
const BASE = {system: 'diploid',             // Inheritance system (changed below)
  neByBlock: false, neBlocks: null,          // Ne is not set separately for each interval
  genTime: 1,                                // Years per generation
  Nf: 1000, Nm: 1000,                        // Breeding females and males per generation
  noDrift: false,                            // Drift is simulated (switched off in check a)
  p0: 0.1,                                   // Starting frequency of the favoured allele
  s: 0,                                      // Selection coefficient (changed below)
  h: 0.5,                                    // Dominance
  onset: 1935,                               // Year selection starts
  years: [1935, 1955, 1975, 2025],           // Sampling years
  n: [20, 20, 20, 40],                       // Individuals sampled at each date
  perfect: true,                             // Perfect genotypes
  depthEarly: 10, depthLast: 15, err: 0.005, minDP: 3,   // Sequencing settings (not used with perfect genotypes)
  alpha: 0.05,                               // Significance level
  m: 1,                                      // One candidate locus
  sMale: 1, decline: 1, mig: 0, pMig: -1,    // No sex difference in selection, constant Ne, no migration
  kSib: 1, maleFrac: 0,                      // Unrelated individuals, females only
  reps: 5000,                                // Selection simulations per run
  seed: 7};                                  // Random seed
const out = {power: [], drift: [], selection: []};   // Results of the three checks
const t0 = Date.now();                       // Starting time, to print the time used


///////////////////////////////////////////////////////////////
// (a) Power without drift and with perfect genotypes        //
// Only sampling is random, so power can be computed exactly //
// (theory_check.py)                                         //
///////////////////////////////////////////////////////////////
// (a) power with no drift and perfect genotypes: only sampling noise remains
const SYSTEMS = [['haploid', 0.5], ['diploid', 0.5], ['diploid', 0], ['diploid', 1], ['haplodiploid', 0.5]];   // [inheritance system, dominance]
const DESIGNS = [['4 points', [1935, 1955, 1975, 2025], [20, 20, 20, 40]],
  ['2 points', [1935, 2025], [20, 40]]];     // [label, sampling years, individuals per date]
let k = 0;                                   // Counter, to give each run its own seed
for (const [system, h] of SYSTEMS)           // Loop over the inheritance systems and dominance values
  for (const p0 of [0.05, 0.2])              // Loop over the starting frequencies
    for (const s of [0.005, 0.01, 0.02, 0.05])   // Loop over the selection coefficients
      for (const [dl, years, n] of DESIGNS) {    // Loop over the two sampling designs
  const P = Object.assign({}, BASE, {system, h, p0, s, years, n, noDrift: true, seed: 100 + k++});   // Settings of this run (drift switched off)
  const r = S.runPower(P);                   // One SPOT run
  out.power.push({system, h, p0, s, design: dl, years, n, reps: P.reps,
    trend: {crit: r.trend.crit, gamma: r.trend.gamma, ties: r.trend.ties, powerEmp: r.trend.powerEmp, powerNaive: r.trend.powerNaive, fprNaive: r.trend.fprNaive},   // Trend test: threshold, ties, powers, false-positive rate
    two: {crit: r.two.crit, gamma: r.two.gamma, ties: r.two.ties, powerEmp: r.two.powerEmp, powerNaive: r.two.powerNaive, fprNaive: r.two.fprNaive}});               // First-versus-last test: the same values
}
console.error('power', ((Date.now() - t0) / 1000).toFixed(0) + 's');   // Progress message with the time used


////////////////////////////////////////////////////////////////////////
// (b) Loss of neutral genetic diversity by drift                     //
// Mean of p(1 - p) relative to p0(1 - p0), at every generation, from //
// 100,000 neutral trajectories                                       //
////////////////////////////////////////////////////////////////////////
// (b) neutral drift: mean of p(1 - p) relative to p0(1 - p0), every generation, from 100,000 neutral trajectories
// overall allele frequency: haploid p; diploid (pf + pm) / 2; haplodiploid (2 pf + pm) / 3 (copies carried by each sex)
const DRIFT = [                              // [label, settings that differ from BASE]
  ['haploid, N = 100', {system: 'haploid', Nf: 100, Nm: 100}],
  ['haploid, N = 1,000', {system: 'haploid', Nf: 1000, Nm: 1000}],
  ['diploid, Nf = Nm = 50', {system: 'diploid', Nf: 50, Nm: 50}],
  ['diploid, Nf = Nm = 500', {system: 'diploid', Nf: 500, Nm: 500}],
  ['diploid, Nf = 100, Nm = 20', {system: 'diploid', Nf: 100, Nm: 20}],
  ['haplodiploid, Nf = Nm = 100', {system: 'haplodiploid', Nf: 100, Nm: 100}],
  ['haplodiploid, Nf = 300, Nm = 50', {system: 'haplodiploid', Nf: 300, Nm: 50}],
  ['diploid, Nf = Nm = 500, decline to 10%', {system: 'diploid', Nf: 500, Nm: 500, decline: 0.1}],
  ['diploid, Ne per interval 1,000 / 50 / 1,000', {system: 'diploid', Nf: 500, Nm: 500, neByBlock: true, neBlocks: [1000, 50, 1000]}],
];
for (const [lab, o] of DRIFT) {              // Loop over the scenarios
  const P = S.prep(Object.assign({}, BASE, o, {p0: 0.3, n: [0, 0, 0, 0], perfect: true}));   // Settings: p0 = 0.3, no individual sampled
  const R = S.makeRng(31),                   // Random number generator
    reps = 100000,                           // Number of neutral trajectories
    G = S.gensOf(P).slice(-1)[0];            // Last generation (90)
  const het = new Float64Array(G + 1),       // Sum of p(1 - p) at each generation
    het2 = new Float64Array(G + 1);          // Sum of its square (for the standard error)
  for (let r = 0; r < reps; r++) {
    const tr = S.simulateRep(P, 0, R, true); // One neutral trajectory (allele frequency in females and in males)
    for (let g = 0; g <= G; g++) {
      const pf = tr.traj[g],                 // Frequency in females at generation g
        pm = tr.trajM[g];                    // Frequency in males at generation g
      const p = o.system === 'haploid' ? pf : o.system === 'diploid' ? (pf + pm) / 2 : (2 * pf + pm) / 3;   // Frequency in the population
      het[g] += p * (1 - p);
      het2[g] += (p * (1 - p)) ** 2;
    }
  }
  out.drift.push({label: lab, params: o, p0: P.p0, G, reps,
    spotHetRetained: S.hetRetained(Object.assign({}, BASE, o)),                   // Value shown by the application
    het: Array.from(het, v => v / reps / (P.p0 * (1 - P.p0))),                    // Simulated diversity left at each generation
    se: Array.from(het, (v, g) => Math.sqrt((het2[g] / reps - (v / reps) ** 2) / reps) / (P.p0 * (1 - P.p0)))});   // Its standard error
}
console.error('drift', ((Date.now() - t0) / 1000).toFixed(0) + 's');   // Progress message with the time used


//////////////////////////////////////////////////////////////
// (c) Trajectories under selection without drift, over 300 //
// generations                                              //
//////////////////////////////////////////////////////////////
// (c) deterministic selection trajectories (no drift), 300 generations
for (const [system, h] of [['haploid', 0.5], ['diploid', 0.5], ['diploid', 0], ['diploid', 1], ['haplodiploid', 0.5], ['haplodiploid', 0]])   // Inheritance systems and dominance values
  for (const s of [0.01, 0.05, 0.2])         // Selection coefficients
    for (const p0 of [0.01, 0.2]) {          // Starting frequencies
    const P = Object.assign({}, BASE, {system, h, s, p0, years: [0, 300], n: [0, 0], onset: 0});   // Settings: 300 generations, no sampling
    const d = S.deterministic(P, true);      // Trajectory without drift (frequency in females and in males)
    out.selection.push({system, h, s, p0, f: Array.from(d.f), m: Array.from(d.m)});
  }

// Mean trajectory with drift in a very large population: it should follow the trajectory without drift
// stochastic mean with a large population should follow the deterministic path (diploid, additive)
{
  const P = S.prep(Object.assign({}, BASE, {system: 'diploid', s: 0.05, p0: 0.01, Nf: 1e6, Nm: 1e6, years: [0, 300], n: [0, 0], onset: 0}));   // Diploid, 1,000,000 breeding females and males
  const R = S.makeRng(5),                    // Random number generator
    reps = 500,                              // Number of trajectories
    mean = new Float64Array(301);            // Mean frequency at each generation
  for (let r = 0; r < reps; r++) {
    const tr = S.simulateRep(P, P.s, R, true);   // One trajectory with selection and drift
    for (let g = 0; g <= 300; g++) mean[g] += tr.traj[g] / reps;
  }
  out.stochasticMean = {system: 'diploid', h: 0.5, s: 0.05, p0: 0.01, N: 1e6, reps, f: Array.from(mean)};
}


//////////////////////////
// Save all the results //
//////////////////////////
fs.writeFileSync('paper/validation/theory_spot.json', JSON.stringify(out));
console.error('done', ((Date.now() - t0) / 1000).toFixed(0) + 's');
