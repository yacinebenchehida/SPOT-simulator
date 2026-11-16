///////////////////////////////////////////////////////////////////////////////////////
// Threshold from a fitted normal distribution, beyond 1,000 candidate loci (S2)     //
//                                                                                   //
// With more than 1,000 candidate loci at alpha = 0.05, SPOT cannot run enough       //
// neutral simulations to read the threshold directly. It then fits a normal         //
// distribution to 200,000 neutral simulations and takes the threshold from it.      //
//                                                                                   //
// This script measures the realised false-positive rate of that threshold, for a    //
// rare (p0 = 0.01) and a common (p0 = 0.1) allele, in the three inheritance         //
// systems (Ne = 1,500; design of the reference case), on 2,000,000 independent      //
// neutral simulations. It also compares the power with the power given by the       //
// threshold read directly from those 2,000,000 simulations (empirical threshold).   //
//                                                                                   //
// Run from the repository root:                                                     //
//   node paper/validation/normal_fit_rare.js > paper/validation/normal_fit_rare.tsv //
///////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////
// Load the simulation code //
//////////////////////////////
const {S, with: W} = require('./base.js');   // S: simulation code; W: reference case with some settings changed


/////////////////////
// Values explored //
/////////////////////
const SYS = [['haploid', 1500, 1500], ['diploid', 750, 750], ['haplodiploid', 1000, 1000]];   // [system, breeding females, breeding males]: Ne = 1,500 in the three systems
const NT = 2000000,                          // Independent neutral simulations used to test each threshold
  MS = [2000, 10000];                        // Numbers of candidate loci


////////////////////////////////
// Header of the results file //
////////////////////////////////
console.log(['system','p0','m','alpha_per_locus','test','threshold_normal_fit','threshold_empirical_2e6','fpr_normal_fit',
  'ratio_to_alpha','power_normal_fit','power_empirical','n_neutral_fit','n_test'].join('\t'));


////////////////////////////////////////////////////////
// For each inheritance system and starting frequency //
////////////////////////////////////////////////////////
let seed = 3000;                             // Each random number generator gets its own seed (3001, 3002, ...)
for (const [system, Nf, Nm] of SYS)          // Loop over the inheritance systems
  for (const p0 of [0.01, 0.1]) {            // Loop over the starting frequencies (rare, common)
  const base = W({system, Nf, Nm, p0, s: 0.05});   // Settings of this case
  const prep = S.prep(S.rescale(base).P),    // Settings in the form used by the simulation code
    t = S.gensOf(prep),                      // Generation number of each sampling date
    K = t.length,                            // Number of sampling dates
    ALL = t.map((_, i) => i),                // Indices of all dates (trend test)
    TWO = [0, K - 1];                        // Indices of the first and last dates (first-versus-last test)

  // Function: score Z computed on the sampling dates listed in idx
  const Z = (x, n, idx) => {
    let N=0,X=0,T=0;                         // Gene copies sampled, copies of the allele, weighted sum of the generation numbers
    for (const i of idx){
      N+=n[i];
      X+=x[i];
      T+=n[i]*t[i]
    }
    if(!N) return 0;                         // No gene copy sampled: Z = 0
    const pb=X/N,                            // Overall allele frequency in the samples
      tb=T/N;                                // Mean generation number of the sampled gene copies
    if(pb<=0||pb>=1) return 0;               // Allele absent or fixed in all samples: Z = 0
    let a=0,b=0;                             // Numerator and sum of squared deviations of time
    for (const i of idx){
      a+=x[i]*(t[i]-tb);
      b+=n[i]*(t[i]-tb)**2
    }
    return b>0 ? a/Math.sqrt(pb*(1-pb)*b) : 0;   // Cochran-Armitage trend score
  };

  // 2,000,000 independent neutral simulations: score Z of both tests
  const R = S.makeRng(++seed),               // Random number generator
    test = {trend: new Float64Array(NT), two: new Float64Array(NT)};   // Z of the trend test and of the first-versus-last test
  for (let i = 0; i < NT; i++) {
    const o = S.simulateRep(prep, 0, R, false);   // One neutral simulation
    test.trend[i] = Z(o.x, o.n, ALL);
    test.two[i] = Z(o.x, o.n, TWO);
  }
  const sorted = {trend: Float64Array.from(test.trend).sort(), two: Float64Array.from(test.two).sort()};   // Sorted copies, to read empirical thresholds

  // For each number of candidate loci: threshold of SPOT (normal fit) against the empirical threshold
  for (const m of MS) {
    const r = S.runPower(Object.assign({}, base, {m, seed: ++seed})),   // One SPOT run with m candidate loci
      a = 0.05 / m;                          // Significance level per locus
    for (const k of ['trend', 'two']) {      // For both tests
      const c = r[k].crit,                   // Threshold of SPOT (normal fit)
        E = S.empiricalThreshold(sorted[k], a),
        ce = E.crit;   // empirical threshold of the 2e6 simulations, ties as in SPOT
      let f = 0;
      for (const z of test[k]) f += S.exceed(z, c, r[k].gamma);         // Neutral simulations above the threshold of SPOT
      let pe = 0;
      for (const z of r.sel[k]) pe += S.exceed(z, ce, E.gamma);         // Selection simulations above the empirical threshold
      console.log([system, p0, m, a.toExponential(1), k, c.toFixed(3), ce.toFixed(3), (f / NT).toExponential(2), ((f / NT) / a).toFixed(1),
        r[k].powerEmp.toFixed(3), (pe / r.sel[k].length).toFixed(3), r.nNeutral, NT].join('\t'));   // One results line
    }
  }
  console.error(system, p0, 'done');         // Progress message on the screen
}
