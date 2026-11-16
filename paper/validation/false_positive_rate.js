/////////////////////////////////////////////////////////////////////////////////////////
// False-positive rate of the two thresholds (Figure 2, Tables S1 to S3)               //
//                                                                                     //
// Question: when there is no selection, how often does each threshold call a          //
// neutral locus "selected"? The intended rate is alpha' = alpha / number of loci.     //
//                                                                                     //
// For each design, SPOT sets its drift-aware threshold from its own neutral           //
// simulations (20 independent runs, with the default 1,000 selection and 2,000        //
// neutral simulations each). The realised false-positive rate of each threshold is    //
// then measured on a large, independent set of neutral simulations. The               //
// sampling-only threshold (the standard normal quantile for alpha') is evaluated      //
// on the same independent set.                                                        //
//                                                                                     //
// Four parts, chosen by a letter on the command line:                                 //
//   A : effective population size x inheritance system         (Figure 2a, Table S1)  //
//   B : one design factor changed at a time                    (Figure 2b, Table S2)  //
//   C : many candidate loci (small alpha per locus)            (Figure 2c, Table S3)  //
//   G : as C, but with the threshold from a normal distribution fitted to only        //
//       2,000 neutral simulations                              (Figure 2c, dashed)    //
//                                                                                     //
// Run from the repository root:                                                       //
//   node paper/validation/false_positive_rate.js A > paper/validation/fpr_A.tsv       //
// Part C can be split to run on two processor cores, for example "C 0,1" and "C 2,3". //
/////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////
// Load the simulation code //
//////////////////////////////
const S = require('../../sim.js');


////////////////////////////////////////////////////////////
// Settings of the reference design (no selection: s = 0) //
////////////////////////////////////////////////////////////
const BASE = {system: 'haplodiploid',        // Inheritance system
  neByBlock: false, neBlocks: null,          // Ne is not set separately for each interval
  genTime: 1,                                // Years per generation
  Nf: 1000, Nm: 1000,                        // Breeding females and males per generation (Ne = 1,500)
  noDrift: false,                            // Drift is simulated
  p0: 0.1,                                   // Starting frequency of the allele
  s: 0,                                      // No selection
  h: 0.5,                                    // Dominance
  onset: 1935,                               // Year selection starts (not used, s = 0)
  years: [1935, 1955, 1975, 2025],           // Sampling years
  n: [20, 20, 20, 40],                       // Individuals sampled at each date
  perfect: false,                            // Genotypes are called from sequencing reads
  depthEarly: 10, depthLast: 15,             // Mean sequencing depth (earlier dates, last date)
  err: 0.005,                                // Error rate per read
  minDP: 3,                                  // Minimum number of reads to call a genotype
  alpha: 0.05,                               // Significance level
  m: 1,                                      // Number of candidate loci
  sMale: 1, decline: 1, mig: 0, pMig: -1,    // No sex difference in selection, constant Ne, no migration
  kSib: 1, maleFrac: 0,                      // Unrelated individuals, females only
  reps: 1000,                                // Selection simulations per run (default)
  seed: 1};                                  // Random seed
const N_THRESHOLDS = 20;                     // Number of independent SPOT runs (one threshold each)


///////////////////////////////////////////////////////////////////////
// Function: the settings P with the inheritance system "system" and //
// the breeding numbers that give the effective population size Ne   //
///////////////////////////////////////////////////////////////////////
// breeding numbers that give a chosen Ne in each system
function withNe(P, system, Ne) {
  const Q = Object.assign({}, P, {system});  // Copy of the settings with the chosen system
  if (system === 'haploid') {                // Haploid: Ne = N
    Q.Nf = Ne;
    Q.Nm = Ne;
  }
  else if (system === 'diploid') {           // Diploid with as many females as males: Ne = Nf + Nm
    Q.Nf = Ne / 2;
    Q.Nm = Ne / 2;
  }
  else {
    Q.Nf = Ne / 1.5;
    Q.Nm = Ne / 1.5;
  }   // haplodiploid with Nf = Nm: Ne = 9N^2 / 6N = 1.5 N
  return Q;
}


///////////////////////////////////////////////////////////////////
// Function: score Z of nWanted independent neutral simulations, //
// for the trend test and for the first-versus-last test         //
///////////////////////////////////////////////////////////////////
// independent neutral Z values: exactly nWanted neutral simulations (s = 0), drawn with SPOT's own simulation code
function neutralZ(P, nWanted, seed) {
  const sc = S.rescale(Object.assign({}, P, {s: 0})),   // Settings without selection (rescaled if the series is very long)
    Q = S.prep(sc.P),                        // The same settings in the form used by the simulation code
    R = S.makeRng(seed);                     // Random number generator
  const t = S.gensOf(sc.P),                  // Generation number of each sampling date
    K = t.length,                            // Number of sampling dates
    ALL = t.map((_, i) => i),                // Indices of all dates (trend test)
    TWO = [0, K - 1];                        // Indices of the first and last dates (first-versus-last test)
  const out = {trend: new Float64Array(nWanted), two: new Float64Array(nWanted)};   // Z of both tests
  for (let i = 0; i < nWanted; i++) {
    const o = S.simulateRep(Q, 0, R, false); // One neutral simulation
    out.trend[i] = trendZ(o.x, o.n, t, ALL);
    out.two[i] = trendZ(o.x, o.n, t, TWO);
  }
  return out;
}


////////////////////////////////////////////////////////////////////////
// Function: score Z of the trend test                                //
// x: copies of the allele per date; n: gene copies sampled per date; //
// t: generation number of each date; idx: dates used by the test     //
////////////////////////////////////////////////////////////////////////
function trendZ(x, n, t, idx) {   // Cochran-Armitage trend statistic, as trendZ() in sim.js
  let N = 0,                                 // Total number of gene copies sampled
    X = 0,                                   // Total number of copies of the allele
    T = 0;                                   // Sum of the generation numbers, weighted by the gene copies sampled
  for (const i of idx) {                     // Loop over the sampling dates
    N += n[i];
    X += x[i];
    T += n[i] * t[i];
  }
  if (N === 0) return 0;                     // No gene copy sampled: Z = 0
  const pb = X / N,                          // Overall allele frequency in the samples
    tb = T / N;                              // Mean generation number of the sampled gene copies
  if (pb <= 0 || pb >= 1) return 0;          // Allele absent or fixed in all samples: Z = 0
  let num = 0,                               // Numerator: covariance between allele copies and time
    den = 0;                                 // Sum of squared deviations of time
  for (const i of idx) {                     // Loop over the sampling dates
    num += x[i] * (t[i] - tb);
    den += n[i] * (t[i] - tb) ** 2;
  }
  return den <= 0 ? 0 : num / Math.sqrt(pb * (1 - pb) * den);
}


////////////////////////////////////////////////////////////////////
// Function: false-positive rates of one design                   //
// Returns two results lines (trend test, first-versus-last test) //
////////////////////////////////////////////////////////////////////
function evaluate(label, P, nTest) {
  // 20 SPOT runs: one drift-aware threshold each
  const thr = [];
  for (let i = 0; i < N_THRESHOLDS; i++)
    thr.push(S.runPower(Object.assign({}, P, {seed: 1000 + i})));

  // Independent neutral simulations, and the sampling-only threshold
  const test = neutralZ(P, nTest, 99991);    // Z of nTest independent neutral simulations
  const a = P.alpha / P.m,                   // Significance level per locus (alpha')
    zc = S.qnorm(1 - a);                     // Sampling-only threshold

  const rows = [];
  for (const k of ['trend', 'two']) {        // For both tests
    const z = test[k],                       // Z of the independent neutral simulations
      nT = z.length;                         // Their number

    // False-positive rate of each of the 20 drift-aware thresholds
    const fpr = thr.map(r => {
      let c = 0;
      for (const v of z) c += S.exceed(v, r[k].crit, r[k].gamma);
      return c / nT;
    });   // ties at the threshold count with weight gamma (sim.js)

    // False-positive rate of the sampling-only threshold
    let naive = 0;
    for (const v of z) if (v > zc) naive++;

    // Summary over the 20 thresholds
    const mean = fpr.reduce((x, y) => x + y, 0) / fpr.length;                                       // Mean false-positive rate
    const sd = Math.sqrt(fpr.reduce((x, y) => x + (y - mean) * (y - mean), 0) / (fpr.length - 1));  // Its standard deviation between runs
    const zero = z.reduce((c, v) => c + (v === 0 ? 1 : 0), 0) / nT;                                 // Proportion of neutral simulations with Z = 0 (allele lost or fixed)

    // One results line
    rows.push([label, k, P.system, Math.round(S.neFor(P)), P.p0, P.years.length, S.gensOf(P).slice(-1)[0], +a.toPrecision(3),
      thr[0][k].method, mean.toExponential(4), sd.toExponential(3), Math.min(...fpr).toExponential(3), Math.max(...fpr).toExponential(3),
      (naive / nT).toExponential(4), nT, (thr.reduce((x, r) => x + r[k].crit, 0) / thr.length).toFixed(3), zc.toFixed(3), zero.toFixed(4)]);
  }
  return rows;
}


/////////////////////////////////////////////////////////////////
// Part chosen on the command line, header and output function //
/////////////////////////////////////////////////////////////////
const part = process.argv[2] || 'A';         // Letter of the part to run (A by default)
const header = ['label', 'statistic', 'system', 'Ne', 'p0', 'points', 'generations', 'alpha_per_locus', 'threshold_method',
  'fpr_drift_aware_mean', 'fpr_drift_aware_sd', 'fpr_drift_aware_min', 'fpr_drift_aware_max', 'fpr_naive', 'n_test',
  'threshold_mean', 'threshold_naive', 'share_Z_zero'];
console.log(header.join('\t'));
const t0 = Date.now();                       // Starting time, to print the time used
const emit = rows => {                       // Function: write the results lines of one design and a progress message
  for (const r of rows) console.log(r.join('\t'));
  console.error(rows[0][0], ((Date.now() - t0) / 1000).toFixed(0) + 's');
};


//////////////////////////////////////////////////////////////////
// Part A: 8 effective population sizes x 3 inheritance systems //
// 100,000 independent neutral simulations per design           //
//////////////////////////////////////////////////////////////////
if (part === 'A') { // effective size x inheritance system
  for (const system of ['haploid', 'diploid', 'haplodiploid'])
    for (const Ne of [100, 300, 1000, 3000, 10000, 30000, 100000, 1000000])
      emit(evaluate(`A ${system} Ne=${Ne}`, withNe(BASE, system, Ne), 100000));
}


////////////////////////////////////////////////////////////////////
// Part B: the reference design with one factor changed at a time //
// [label, settings that differ from the reference design]        //
////////////////////////////////////////////////////////////////////
if (part === 'B') { // one design factor at a time, haplodiploid, Ne = 1,500 (1,000 breeding females and males)
  const ten = [1935, 1945, 1955, 1965, 1975, 1985, 1995, 2005, 2015, 2025];   // Ten sampling years, one every ten years
  const V = [
    ['base: 4 points, 20/20/20/40, 10x/15x, p0 = 0.1', {}],
    ['2 points (1935, 2025)', {years: [1935, 2025], n: [20, 40]}],
    ['10 points, 20 each', {years: ten, n: ten.map(() => 20)}],
    ['p0 = 0.01', {p0: 0.01}],
    ['p0 = 0.5', {p0: 0.5}],
    ['depth 2x at all points', {depthEarly: 2, depthLast: 2, minDP: 1}],
    ['perfect genotypes', {perfect: true}],
    ['siblings: 4 per family', {kSib: 4}],
    ['decline to 10% by the last point', {decline: 0.1}],
    ['bottleneck: Ne 1,500 / 100 / 1,500', {neByBlock: true, neBlocks: [1500, 100, 1500]}],
    ['migration 1% per generation', {mig: 0.01}],
    ['200 individuals per point', {n: [200, 200, 200, 200]}],
  ];
  for (const [lab, o] of V)
    emit(evaluate('B ' + lab, Object.assign({}, BASE, o), 100000));
}


//////////////////////////////////////////////////////////////////////
// Part C: 1, 10, 100 and 1,000 candidate loci, for four designs    //
// SPOT runs more neutral simulations when alpha per locus is small //
//////////////////////////////////////////////////////////////////////
if (part === 'C') { // small per-locus alpha (Bonferroni): SPOT raises the number of neutral replicates (neutralReps in sim.js)
  const D = [['haplodiploid Ne=300', withNe(BASE, 'haplodiploid', 300)],
    ['haplodiploid Ne=1500', BASE],
    ['haplodiploid Ne=100000', withNe(BASE, 'haplodiploid', 100000)],
    ['haplodiploid Ne=1500, p0=0.01', Object.assign({}, BASE, {p0: 0.01})]];   // The four designs
  const pick = process.argv[3] ? process.argv[3].split(',').map(Number) : D.map((_, i) => i);   // optional subset of designs, e.g. 0,1
  for (const [lab, P] of pick.map(i => D[i]))
    for (const m of [1, 10, 100, 1000])                                        // Numbers of candidate loci
      emit(evaluate(`C ${lab} m=${m}`, Object.assign({}, P, {m}), m >= 100 ? 1000000 : 100000));   // 1,000,000 independent neutral simulations from 100 loci
}


////////////////////////////////////////////////////////////////////////
// Part G: as part C, but SPOT may use only 2,000 neutral simulations //
// The threshold then comes from a fitted normal distribution         //
////////////////////////////////////////////////////////////////////////
if (part === 'G') { // comparison: threshold from a normal distribution fitted to 2,000 neutral simulations (neutralMax = 2,000),
                    // i.e. without the extra neutral simulations of part C (Gaussian fit as soon as alpha' x 2,000 < 10)
  const D = [['haplodiploid Ne=300', withNe(BASE, 'haplodiploid', 300)],
    ['haplodiploid Ne=1500', BASE],
    ['haplodiploid Ne=100000', withNe(BASE, 'haplodiploid', 100000)],
    ['haplodiploid Ne=1500, p0=0.01', Object.assign({}, BASE, {p0: 0.01})]];   // The same four designs
  const pick = process.argv[3] ? process.argv[3].split(',').map(Number) : D.map((_, i) => i);
  for (const [lab, P] of pick.map(i => D[i]))
    for (const m of [1, 10, 100, 1000])
      emit(evaluate(`G ${lab} m=${m}`, Object.assign({}, P, {m, neutralMax: 2000}), m >= 100 ? 1000000 : 100000));
}
