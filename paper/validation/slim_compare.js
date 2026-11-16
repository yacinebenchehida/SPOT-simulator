///////////////////////////////////////////////////////////////////////////////////////
// Comparison of SPOT with SLiM: allele frequencies, power, false positives          //
// (Supporting Information S4, Table S5)                                             //
//                                                                                   //
// For each scenario of slim_runs.sh, SPOT simulates the same number of populations  //
// (2,000) with the same settings, with its own code (simulateRep in sim.js). Both   //
// programs record, at generations 0, 20, 40 and 90:                                 //
//   pfs : allele frequency among breeding females after selection (haploid: all     //
//         individuals)                                                              //
//   x/n : copies of the allele among the sampled females (perfect genotypes)        //
//                                                                                   //
// Then, separately for each program and with exactly the same code:                 //
//   1. the distribution of pfs is summarised (mean, SD) and the two distributions   //
//      are compared (difference of means in standard errors, ratio of SDs, and a    //
//      Kolmogorov-Smirnov test at generation 90);                                   //
//   2. the score Z of the trend test is computed for every simulation; the          //
//      drift-aware threshold is the 95% quantile of Z in the neutral simulations    //
//      (as in SPOT), and power is the proportion of selection simulations above it; //
//   3. the threshold of SPOT is also applied to the neutral simulations of SLiM:    //
//      if SPOT simulates drift correctly, this false-positive rate should be close  //
//      to 5%. For a precise version of this check (slim_fpr.tsv), the threshold of  //
//      SPOT is taken from 20,000 neutral simulations and applied to all the neutral //
//      simulations of SLiM for the scenario, including the extra runs _b, _c, _d    //
//      of slim_runs.sh when present (up to 8,000).                                  //
//                                                                                   //
// Run from the repository root, after paper/validation/slim_runs.sh:                //
//   node paper/validation/slim_compare.js                                           //
// It writes paper/validation/slim_results.tsv (frequencies), slim_power.tsv (power) //
// and slim_fpr.tsv (false positives).                                               //
///////////////////////////////////////////////////////////////////////////////////////


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
  R = 2000;                                  // Simulations per scenario


/////////////////////////////////////////////////////////////////////////
// The six scenarios with selection                                    //
// [name of the SLiM run with selection, name of the neutral SLiM run, //
//  system, breeding females, breeding males, p0, s, dominance]        //
/////////////////////////////////////////////////////////////////////////
// same scenarios as slim_runs.sh: [selected run, matching neutral run, system, NF, NM, P0, S, H]
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
const over = (z, c) => S.exceed(z, c.crit, c.gamma);   // 1 above, gamma when tied, 0 below

const mean = a => a.reduce((p, q) => p + q, 0) / a.length;   // Mean of a list of numbers

const sd = a => {                            // Standard deviation of a list of numbers
  const m = mean(a);
  return Math.sqrt(a.reduce((p, q) => p + (q - m) ** 2, 0) / (a.length - 1));
};

// Kolmogorov-Smirnov test for two samples: are the two distributions the same?
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


///////////////////////////////////////////////////////////////////////
// Parts 1 and 2: allele frequencies and power, scenario by scenario //
///////////////////////////////////////////////////////////////////////
// ---- comparison ----
const rows = [['scenario', 'run', 'generation', 'mean_spot', 'mean_slim', 'z_mean', 'sd_spot', 'sd_slim', 'sd_ratio', 'ks_D', 'ks_p']];   // Header of the table of frequencies
const prow = [['scenario', 'system', 'NF', 'NM', 'p0', 's', 'h', 'threshold_spot', 'threshold_slim', 'power_spot', 'power_slim',
  'se_difference', 'z_power', 'fpr_spot_threshold_on_slim_neutral']];   // Header of the table of power
const done = new Set(),                      // Runs whose frequencies are already compared (a neutral run can serve two scenarios)
  neutralSpec = new Map();                   // Settings of each neutral scenario, kept for part 3
SC.forEach(([sel, neu, system, Nf, Nm, p0, s, h], k) => {
  neutralSpec.set(neu, [system, Nf, Nm, p0, h]);

  // Simulations of both programs: [SLiM, SPOT], with selection and neutral
  const runs = {sel: [readSlim(sel), runSpot(system, Nf, Nm, p0, s, h, 100 + k)],
                neu: [readSlim(neu), runSpot(system, Nf, Nm, p0, 0, h, 200 + k)]};

  // Part 1: allele frequency of the population at generations 20, 40 and 90
  for (const [lab, name] of [['selected', sel], ['neutral', neu]]) {
    if (done.has(name)) continue;            // Skip a run that is already compared
    done.add(name);
    const [sl, sp] = runs[lab === 'selected' ? 'sel' : 'neu'];   // SLiM and SPOT simulations of this run
    GENS.forEach((g, i) => {
      if (g === 0) return;                   // Generation 0: the frequency is p0 in both programs
      const a = sp.map(r => r.pfs[i]),       // Frequencies in SPOT
        b = sl.map(r => r.pfs[i]);           // Frequencies in SLiM
      const z = (mean(a) - mean(b)) / Math.sqrt(sd(a) ** 2 / a.length + sd(b) ** 2 / b.length);   // Difference of the means, in standard errors
      const [D, p] = g === 90 ? ks(a, b) : ['', ''];   // Kolmogorov-Smirnov test at generation 90 only
      rows.push([name, lab, g, mean(a).toFixed(4), mean(b).toFixed(4), z.toFixed(2), sd(a).toFixed(4), sd(b).toFixed(4),
        (sd(a) / sd(b)).toFixed(3), D === '' ? '' : D.toFixed(4), p === '' ? '' : p.toFixed(3)]);
    });
  }

  // Part 2: power, each program with its own threshold
  const Z = rs => rs.map(r => trendZ(r.x, r.n));   // Score Z of each simulation of a list
  const [zSelSlim, zSelSpot, zNeuSlim, zNeuSpot] = [Z(runs.sel[0]), Z(runs.sel[1]), Z(runs.neu[0]), Z(runs.neu[1])];
  const cSpot = quant95(zNeuSpot),           // Threshold from the neutral simulations of SPOT
    cSlim = quant95(zNeuSlim);               // Threshold from the neutral simulations of SLiM
  const pSpot = mean(zSelSpot.map(z => over(z, cSpot))),   // Power in SPOT
    pSlim = mean(zSelSlim.map(z => over(z, cSlim)));       // Power in SLiM
  const se = Math.sqrt(pSpot * (1 - pSpot) / R + pSlim * (1 - pSlim) / R);   // Standard error of the difference between the two powers
  const fpr = mean(zNeuSlim.map(z => over(z, cSpot)));     // False-positive rate of the threshold of SPOT on the neutral simulations of SLiM
  prow.push([sel, system, Nf, Nm, p0, s, h, cSpot.crit.toFixed(3), cSlim.crit.toFixed(3), pSpot.toFixed(4), pSlim.toFixed(4),
    se.toFixed(4), ((pSpot - pSlim) / se).toFixed(2), fpr.toFixed(4)]);
  console.error(sel, 'power SPOT', pSpot.toFixed(3), 'SLiM', pSlim.toFixed(3), 'FPR of SPOT threshold on SLiM neutral', fpr.toFixed(3));   // Progress message
});
fs.writeFileSync('paper/validation/slim_results.tsv', rows.map(r => r.join('\t')).join('\n') + '\n');   // Write the table of frequencies
fs.writeFileSync('paper/validation/slim_power.tsv', prow.map(r => r.join('\t')).join('\n') + '\n');     // Write the table of power


/////////////////////////////////////////////////////////////////
// Part 3: false-positive rate of the threshold of SPOT on the //
// neutral populations of SLiM (precise version)               //
/////////////////////////////////////////////////////////////////
// ---- precise check of the false-positive rate of SPOT's threshold on SLiM's neutral populations ----
const frow = [['neutral_scenario', 'system', 'NF', 'NM', 'p0', 'threshold_spot_20000', 'n_slim', 'fpr_on_slim', 'se', 'z_vs_5pc',
  'sd_Z_spot', 'sd_Z_slim', 'fpr_on_independent_spot']];   // Header
let k2 = 0;                                  // Counter, to give each scenario its own seeds
for (const [neu, [system, Nf, Nm, p0, h]] of neutralSpec) {   // Loop over the neutral scenarios
  const files = [neu, neu + '_b', neu + '_c', neu + '_d'].filter(f => fs.existsSync(`paper/validation/slim_out/${f}.tsv`));   // SLiM neutral runs that exist for this scenario
  const zSlim = files.flatMap(f => readSlim(f).map(r => trendZ(r.x, r.n)));   // Z of all the neutral simulations of SLiM

  // Threshold of SPOT from 20,000 neutral simulations (10 x 2,000)
  const big = [];
  for (let j = 0; j < 10; j++) big.push(...runSpot(system, Nf, Nm, p0, 0, h, 500 + 10 * k2 + j));
  const zSpot = big.map(r => trendZ(r.x, r.n)),
    c = quant95(zSpot);

  // The same threshold on 8,000 independent neutral simulations of SPOT (4 x 2,000)
  const test = [];
  for (let j = 0; j < 4; j++) test.push(...runSpot(system, Nf, Nm, p0, 0, h, 900 + 10 * k2 + j));
  const fprSpot = mean(test.map(r => over(trendZ(r.x, r.n), c)));

  // False-positive rate on the neutral simulations of SLiM, and its standard error if the true rate is 5%
  const fpr = mean(zSlim.map(z => over(z, c))),
    se = Math.sqrt(0.05 * 0.95 / zSlim.length);
  frow.push([neu, system, Nf, Nm, p0, c.crit.toFixed(3), zSlim.length, fpr.toFixed(4), se.toFixed(4), ((fpr - 0.05) / se).toFixed(2),
    sd(zSpot).toFixed(3), sd(zSlim).toFixed(3), fprSpot.toFixed(4)]);
  console.error(neu, 'SLiM neutral replicates', zSlim.length, 'FPR of SPOT threshold', fpr.toFixed(4), '(SPOT on itself', fprSpot.toFixed(4) + ')');   // Progress message
  k2++;
}
fs.writeFileSync('paper/validation/slim_fpr.tsv', frow.map(r => r.join('\t')).join('\n') + '\n');   // Write the table of false-positive rates
