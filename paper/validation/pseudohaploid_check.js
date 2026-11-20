//////////////////////////////////////////////////////////////////////////////////////////
// Robustness check: power of the human ancient DNA case study with pseudo-haploid calls //
//                                                                                      //
// Process, for each simulated data set:                                                //
//   1. Population: drift and selection as in SPOT (sim.js).                            //
//   2. Sampling: the individuals of each date are drawn from the population.           //
//   3. Sequencing: each individual gets a Poisson number of reads at the locus         //
//      (mean = depth); some individuals get no read.                                   //
//   4. The same reads are analysed in two ways:                                        //
//      - SPOT: each individual with at least one read gets a diploid genotype          //
//        (maximum likelihood) and counts as two gene copies;                           //
//      - pseudo-haploid: one read is drawn at random for each individual with at       //
//        least one read and counts as one gene copy.                                   //
//      An individual without a read gives nothing in either analysis.                  //
//   5. Score Z of the trend test, for each analysis.                                   //
// For each analysis, the threshold is set from neutral simulations analysed the same   //
// way, and power is the share of selection simulations above it.                       //
//                                                                                      //
// Settings: human ancient DNA case study (as in worked_examples.js), with the ancient  //
// samples at 1x, 2x, 5x or 10x (present-day sample at 20x) and s = 0.005 or 0.0075;     //
// 5,000 selection and 10,000 neutral simulations per design.                           //
//                                                                                      //
// The copy of the simulation is checked first: with the pseudo-haploid analysis       //
// switched off, it gives exactly the same power as sim.js.                             //
//                                                                                      //
// Run from the repository root:                                                        //
//   node paper/validation/pseudohaploid_check.js > paper/validation/pseudohaploid.tsv  //
//////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////
// Load the simulation code //
//////////////////////////////
const S = require('../../sim.js');             // Same code as the one used by the application


//////////////////////////////////////////////////////////////////////////////////////
// Function: one simulation, copied from simulateRep in sim.js (diploid individuals, //
// unrelated, females only, as in the human case study). The reads of each sampled   //
// individual are analysed twice: diploid genotype (as in sim.js) and one read drawn //
// at random (pseudo-haploid). With both = false, the second analysis is skipped,    //
// so that the random numbers are exactly those of sim.js (check below).             //
//////////////////////////////////////////////////////////////////////////////////////
function simRep(P, s, R, both) {
  if (P.system !== 'diploid' || P.kSib > 1 || (P.maleFrac || 0) > 0 || P.mig > 0) throw new Error('this check handles unrelated diploid females without migration only');
  const gens = S.gensOf(P), K = gens.length, G = gens[K - 1];   // Sampling generations and last generation
  const onset = (P.onset - P.years[0]) / S.genTimeOf(P), h = P.h;   // Generation selection starts; dominance
  let pf = P.p0, pm = P.p0, si = 0;            // Allele frequency in breeding females and males; next sampling date
  const xd = new Array(K).fill(0), nd = new Array(K).fill(0);   // SPOT analysis: copies of the allele and gene copies, per date
  const xp = new Array(K).fill(0), np = new Array(K).fill(0);   // Pseudo-haploid analysis: the same
  for (let g = 0; g <= G; g++) {               // Loop over generations, as in sim.js
    const sg = g > onset ? s : 0;              // Selection coefficient in this generation
    const AA = pf * pm, Aa = pf * (1 - pm) + pm * (1 - pf);   // Genotype frequencies in offspring before selection
    const w = AA * (1 + sg) + Aa * (1 + h * sg) + (1 - AA - Aa);   // Mean fitness
    const fAA = AA * (1 + sg) / w, fAa = Aa * (1 + h * sg) / w;   // Genotype frequencies after selection
    const pfs = fAA + fAa / 2, pms = pfs;      // Allele frequency after selection (same in both sexes here)
    while (si < K && g === gens[si]) {         // Sampling at this generation
      const last = si === K - 1, dep = last ? P.depthLast : P.depthEarly;   // Mean depth of this date
      const e = last ? P.errLast : P.err, lE = last ? P.lEL : P.lE, l1E = last ? P.l1EL : P.l1E;   // Error rate and its logarithms
      for (let i = 0; i < P.n[si]; i++) {      // Each sampled individual
        const r = R.u(), gt = r < fAA ? 2 : (r < fAA + fAa ? 1 : 0);   // True genotype (copies of the favoured allele)
        if (P.perfect) { xd[si] += gt; nd[si] += 2; continue; }
        const d = R.pois(dep);                 // Number of reads at the locus
        if (d < P.minDP) continue;             // No read (minimum depth 1): nothing in either analysis
        const q = gt === 0 ? e : gt === 1 ? 0.5 : 1 - e;   // Probability that a read shows the favoured allele
        const k = R.binom(d, q);               // Reads showing the favoured allele
        const l0 = k * lE + (d - k) * l1E, l1 = -d * Math.LN2, l2 = k * l1E + (d - k) * lE;   // Log-likelihoods of the three genotypes
        xd[si] += (l1 >= l0 && l1 >= l2) ? 1 : (l2 > l0 ? 2 : 0); nd[si] += 2;   // SPOT: maximum-likelihood genotype, two gene copies
        if (both) { xp[si] += R.u() < k / d ? 1 : 0; np[si] += 1; }   // Pseudo-haploid: one of the d reads at random, one gene copy
      }
      si++;
    }
    const f = S.sizeFactor(P, gens, g), Nf = Math.max(2, Math.round(P.Nf * f)), Nm = Math.max(2, Math.round(P.Nm * f));   // Drift, as in sim.js
    pf = R.binom(2 * Nf, pfs) / (2 * Nf);
    pm = R.binom(2 * Nm, pms) / (2 * Nm);
  }
  return {xd, nd, xp, np};
}


///////////////////////////////////////////////////////////
// Function: score Z of the trend test (as in sim.js)    //
///////////////////////////////////////////////////////////
function trendZ(x, n, t) {
  let N = 0, X = 0, T = 0;                     // Gene copies, copies of the allele, time weighted by gene copies
  for (let i = 0; i < t.length; i++) { N += n[i]; X += x[i]; T += n[i] * t[i]; }
  if (N === 0) return 0;
  const pb = X / N, tb = T / N;                // Overall frequency and mean time
  if (pb <= 0 || pb >= 1) return 0;
  let num = 0, den = 0;
  for (let i = 0; i < t.length; i++) { num += x[i] * (t[i] - tb); den += n[i] * (t[i] - tb) ** 2; }
  return den <= 0 ? 0 : num / Math.sqrt(pb * (1 - pb) * den);
}


////////////////////////////////////////////////////////////////////////
// Function: power of the trend test for both analyses, as runPower   //
// in sim.js (selection simulations first, then neutral ones)         //
////////////////////////////////////////////////////////////////////////
function power(P0, both) {
  const P = S.prep(S.rescale(P0).P), R = S.makeRng(P.seed), nn = S.neutralReps(P), t = S.gensOf(P);
  const Z = {spot: [new Float64Array(P.reps), new Float64Array(nn)], pseudo: [new Float64Array(P.reps), new Float64Array(nn)]};   // Z of selection and neutral simulations
  for (const [k, reps, s] of [[0, P.reps, P.s], [1, nn, 0]])   // Selection simulations, then neutral simulations
    for (let r = 0; r < reps; r++) {
      const o = simRep(P, s, R, both);
      Z.spot[k][r] = trendZ(o.xd, o.nd, t);
      if (both) Z.pseudo[k][r] = trendZ(o.xp, o.np, t);
    }
  const out = {};
  for (const a of both ? ['spot', 'pseudo'] : ['spot']) {
    const {crit, gamma} = S.empiricalThreshold(Z[a][1], P.alpha / P.m);   // Drift-aware threshold from the neutral simulations
    let pe = 0;
    for (const z of Z[a][0]) pe += S.exceed(z, crit, gamma);
    out[a] = {power: pe / P.reps, crit};
  }
  return out;
}


////////////////////////////////////////////////////////////////////////////
// Human ancient DNA case study (settings of worked_examples.js)          //
////////////////////////////////////////////////////////////////////////////
const HUMAN = {system: 'diploid', neByBlock: false, neBlocks: null, genTime: 28, Nf: 5000, Nm: 5000, noDrift: false, p0: 0.1, s: 0.005, h: 0.5,
  onset: -7000, years: [-7000, -5000, -3000, -1000, 1000, 2000], n: [50, 50, 50, 50, 50, 100], perfect: false, depthEarly: 2, depthLast: 20,
  err: 0.01, errLast: 0.01, minDP: 1, alpha: 0.05, m: 1, sMale: 1, decline: 1, mig: 0, pMig: -1, kSib: 1, maleFrac: 0, reps: 5000, seed: 11};


//////////////////////////////////////////////////////////////////////
// Check of the copy: SPOT analysis alone gives the power of sim.js //
//////////////////////////////////////////////////////////////////////
for (const d of [1, 10]) {
  const P = Object.assign({}, HUMAN, {depthEarly: d, reps: 1000, seed: 5});
  const a = power(P, false).spot.power, b = S.runPower(P).trend.powerEmp;
  if (a !== b) throw new Error(`copy of simulateRep differs from sim.js: ${a} vs ${b}`);
}


///////////////////////////////////////////////////////
// Results: one line per design and analysis         //
///////////////////////////////////////////////////////
console.log(['s', 'depth_ancient', 'depth_present_day', 'analysis', 'power_trend', 'threshold'].join('\t'));   // Header
let seed = 3000;                               // One seed per design (3001, 3002, ...), shared by both analyses
for (const s of [0.005, 0.0075])               // Two selection coefficients
  for (const d of [1, 2, 5, 10]) {             // Depth of the ancient samples
    const r = power(Object.assign({}, HUMAN, {s, depthEarly: d, seed: ++seed}), true);
    for (const [a, name] of [['spot', 'SPOT (diploid genotypes)'], ['pseudo', 'pseudo-haploid']])
      console.log([s, d, HUMAN.depthLast, name, r[a].power.toFixed(4), r[a].crit.toFixed(3)].join('\t'));
  }
