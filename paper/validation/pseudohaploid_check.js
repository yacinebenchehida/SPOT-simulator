////////////////////////////////////////////////////////////////////////////////////////////
// Robustness check: do the depth conclusions hold with pseudo-haploid calls?             //
//                                                                                        //
// SPOT calls a diploid genotype for each sampled individual from its reads (maximum      //
// likelihood). Ancient DNA is often analysed with pseudo-haploid calls instead: one read //
// drawn at random per individual, which then counts as a single gene copy. This script   //
// recomputes power with pseudo-haploid calls for:                                        //
//   A : the depth analysis of the museum insect case study (Figure 3e): 1x to 30x at     //
//       every date, Ne = 1,500 and 100,000, p0 = 0.1, s = 0.02                           //
//   B : the human ancient DNA budget comparison (Table S10): 100 reads per ancient date  //
//       spread over 100, 50, 25, 20 or 10 individuals, s = 0.005 and 0.0075             //
// Each design is run with SPOT's own calling (sim.js, unchanged) and with pseudo-haploid //
// calls (the copy of simulateRep below), with the same seed and 5,000 selection          //
// simulations. The thresholds are set from neutral simulations with the same calling.    //
//                                                                                        //
// The copy of simulateRep is checked first: without pseudo-haploid calls it must give    //
// exactly the same power as sim.js.                                                      //
//                                                                                        //
// Run from the repository root:                                                          //
//   node paper/validation/pseudohaploid_check.js > paper/validation/pseudohaploid.tsv    //
////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////
// Load the simulation code //
//////////////////////////////
const {S, with: W} = require('./base.js');   // S: simulation code; W: reference case with some settings changed


//////////////////////////////////////////////////////////////////////////////////////
// Function: one simulation, copied from simulateRep in sim.js, with one change:    //
// with pseudo = true, each sampled individual with at least one read gives one     //
// read drawn at random, counted as one gene copy (allele shown by that read)       //
//////////////////////////////////////////////////////////////////////////////////////
function simRep(P, s, R, pseudo) {
  const Yr = P.years, gens = S.gensOf(P), K = gens.length, G = gens[K - 1];   // Sampling generations and last generation
  const onset = (P.onset - Yr[0]) / S.genTimeOf(P), hap = P.system === 'haplodiploid', hpl = P.system === 'haploid', h = P.h, sMrel = P.sMale;
  const inf = P.N === Infinity, mig = P.mig, pMig = P.pMig;
  let pf = P.p0, pm = P.p0, si = 0;            // Allele frequency in breeding females and males; next sampling date
  const x = new Array(K).fill(0), n = new Array(K).fill(0);   // Copies of the allele and gene copies, per date
  let missing = 0, sampled = 0;                // Individuals without a genotype; individuals sampled
  const mf = P.maleFrac || 0;                  // Probability that a sampled individual is a male
  const W2 = pseudo ? 1 : 2;                   // Gene copies counted for a called diploid individual
  const callOne = (gt, dep, last) => {         // Diploid individual: called genotype (or allele) or -1 (missing), as in sim.js
    sampled++;
    if (P.perfect) return gt;
    const d = R.pois(dep);                     // Number of reads
    if (d < P.minDP) { missing++; return -1; }
    const e = last ? P.errLast : P.err, lE = last ? P.lEL : P.lE, l1E = last ? P.l1EL : P.l1E;
    const q = gt === 0 ? e : gt === 1 ? 0.5 : 1 - e;   // Probability that a read shows the favoured allele
    if (pseudo) return R.u() < q ? 1 : 0;      // Pseudo-haploid call: the allele of one read drawn at random
    const k = R.binom(d, q);                   // Reads showing the favoured allele
    const l0 = k * lE + (d - k) * l1E, l1 = -d * Math.LN2, l2 = k * l1E + (d - k) * lE;   // Log-likelihoods of the three genotypes
    return (l1 >= l0 && l1 >= l2) ? 1 : (l2 > l0 ? 2 : 0);   // Maximum-likelihood genotype (heterozygote in case of a tie)
  };
  const callHap = (a, dep, last) => {          // Haploid individual or male: called allele or -1, as in sim.js
    sampled++;
    if (P.perfect) return a;
    const d = R.pois(dep);
    if (d < P.minDP) { missing++; return -1; }
    const e = last ? P.errLast : P.err;
    if (pseudo) return R.u() < (a ? 1 - e : e) ? 1 : 0;   // Pseudo-haploid call: one read drawn at random
    const k = R.binom(d, a ? 1 - e : e);
    return 2 * k > d ? 1 : (2 * k < d ? 0 : (R.u() < 0.5 ? 1 : 0));   // Majority of the reads (ties at random)
  };
  for (let g = 0; g <= G; g++) {               // Loop over generations, as in sim.js
    const sg = g > onset ? s : 0, sm = sg * sMrel;
    let pfs, pms, fAA = 0, fAa = 0, mAA = 0, mAa = 0;
    if (hpl) { pfs = pf * (1 + sg) / (1 + pf * sg); pms = pfs; }
    else {
      const AA = pf * pm, Aa = pf * (1 - pm) + pm * (1 - pf);
      const w = AA * (1 + sg) + Aa * (1 + h * sg) + (1 - AA - Aa);
      fAA = AA * (1 + sg) / w; fAa = Aa * (1 + h * sg) / w; pfs = fAA + fAa / 2;
      if (hap) { pms = pf * (1 + sm) / (pf * (1 + sm) + 1 - pf); }
      else { const wm = AA * (1 + sm) + Aa * (1 + h * sm) + (1 - AA - Aa); mAA = AA * (1 + sm) / wm; mAa = Aa * (1 + h * sm) / wm; pms = mAA + mAa / 2; }
    }
    while (si < K && g === gens[si]) {         // Sampling at this generation (unrelated individuals only)
      const ni = P.n[si], last = si === K - 1, dep = last ? P.depthLast : P.depthEarly;
      if (P.kSib > 1) throw new Error('relatives are not handled in this check');
      for (let i = 0; i < ni; i++) {
        if (hpl) { const c = callHap(R.u() < pfs ? 1 : 0, dep, last); if (c >= 0) { x[si] += c; n[si] += 1; } continue; }
        if (mf > 0 && R.u() < mf) {            // Male
          if (hap) { const c = callHap(R.u() < pms ? 1 : 0, dep, last); if (c >= 0) { x[si] += c; n[si] += 1; } }
          else { const r = R.u(); const gt = r < mAA ? 2 : (r < mAA + mAa ? 1 : 0); const c = callOne(gt, dep, last); if (c >= 0) { x[si] += c; n[si] += W2; } }
        } else {                               // Female (or diploid individual)
          const r = R.u(); const gt = r < fAA ? 2 : (r < fAA + fAa ? 1 : 0); const c = callOne(gt, dep, last); if (c >= 0) { x[si] += c; n[si] += W2; }
        }
      }
      si++;
    }
    if (inf) { pf = pfs; pm = pms; }            // Drift, as in sim.js
    else {
      const f = S.sizeFactor(P, gens, g), Nf = Math.max(2, Math.round(P.Nf * f)), Nm = Math.max(2, Math.round(P.Nm * f));
      if (hpl) { pf = R.binom(Nf, pfs) / Nf; pm = pf; }
      else {
        pf = R.binom(2 * Nf, pfs) / (2 * Nf);
        pm = hap ? R.binom(Nm, pms) / Nm : R.binom(2 * Nm, pms) / (2 * Nm);
      }
    }
    if (mig > 0) { const pmg = pMig < 0 ? P.p0 : pMig; pf = (1 - mig) * pf + mig * pmg; pm = (1 - mig) * pm + mig * pmg; }   // Migration
  }
  return {x, n, missing, sampled};
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


///////////////////////////////////////////////////////////////////////////
// Function: power of the trend test, as runPower in sim.js (selection   //
// simulations first, then neutral ones, same random number generator)   //
///////////////////////////////////////////////////////////////////////////
function power(P0, pseudo) {
  const P = S.prep(S.rescale(P0).P), R = S.makeRng(P.seed), nn = S.neutralReps(P), t = S.gensOf(P);
  const sel = new Float64Array(P.reps), neu = new Float64Array(nn);
  let miss = 0, samp = 0;
  for (let r = 0; r < P.reps; r++) { const o = simRep(P, P.s, R, pseudo); sel[r] = trendZ(o.x, o.n, t); miss += o.missing; samp += o.sampled; }
  for (let r = 0; r < nn; r++) { const o = simRep(P, 0, R, pseudo); neu[r] = trendZ(o.x, o.n, t); }
  const a = P.alpha / P.m;
  const {crit, gamma} = S.empiricalThreshold(neu, a);   // Drift-aware threshold from the neutral simulations
  let pe = 0;
  for (const z of sel) pe += S.exceed(z, crit, gamma);
  return {power: pe / P.reps, crit, missing: miss / samp};
}


/////////////////////////////////////////////////////////////////////
// Check of the copy: without pseudo-haploid calls, same power as  //
// sim.js, to the last digit                                       //
/////////////////////////////////////////////////////////////////////
const HUMAN = {system: 'diploid', neByBlock: false, neBlocks: null, genTime: 28, Nf: 5000, Nm: 5000, noDrift: false, p0: 0.1, s: 0.005, h: 0.5,
  onset: -7000, years: [-7000, -5000, -3000, -1000, 1000, 2000], n: [50, 50, 50, 50, 50, 100], perfect: false, depthEarly: 2, depthLast: 20,
  err: 0.01, minDP: 1, alpha: 0.05, m: 1, sMale: 1, decline: 1, mig: 0, pMig: -1, kSib: 1, maleFrac: 0, reps: 5000, seed: 11};   // Human ancient DNA example (worked_examples.js)
for (const P of [W({depthEarly: 2, depthLast: 2, minDP: 1, reps: 1000, seed: 5}), Object.assign({}, HUMAN, {reps: 1000, seed: 5})]) {
  const a = power(P, false).power, b = S.runPower(P).trend.powerEmp;
  if (a !== b) throw new Error(`copy of simulateRep differs from sim.js: ${a} vs ${b}`);
}


////////////////////////////////////////////////////////////
// Results: one line per design, with both kinds of calls //
////////////////////////////////////////////////////////////
console.log(['part', 'design', 'Ne', 's', 'individuals_per_date', 'depth', 'calls', 'power_trend', 'threshold', 'share_missing'].join('\t'));   // Header
const out = (part, design, Ne, s, n, dep, P) => {   // Run one design with both calls and print two lines
  for (const pseudo of [false, true]) {
    const r = power(P, pseudo);
    console.log([part, design, Ne, s, n, dep, pseudo ? 'pseudo-haploid' : 'SPOT (diploid genotypes)', r.power.toFixed(4), r.crit.toFixed(3), r.missing.toFixed(4)].join('\t'));
  }
};
let seed = 2000;                               // Each design gets its own seed (2001, 2002, ...), the same for both calls
for (const [Nf, Ne] of [[1000, 1500], [66667, 100000]])   // A: museum insect, two levels of drift
  for (const d of [1, 2, 5, 10, 30])           // Same depth at every date
    out('A', `${d}x all dates`, Ne, 0.02, '20/20/20/40', d, W({Nf, Nm: Nf, s: 0.02, depthEarly: d, depthLast: d, minDP: 1, reps: 5000, seed: ++seed}));
for (const s of [0.005, 0.0075])               // B: human ancient DNA, 100 reads per ancient date
  for (const [n, d] of [[100, 1], [50, 2], [25, 4], [20, 5], [10, 10]])
    out('B', `${n} individuals at ${d}x`, 10000, s, n, d, Object.assign({}, HUMAN, {n: [n, n, n, n, n, 100], depthEarly: d, s, seed: ++seed}));
