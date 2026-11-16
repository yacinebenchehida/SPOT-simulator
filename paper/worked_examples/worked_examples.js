////////////////////////////////////////////////////////////////////////////////////////////
// Case studies of the article: what SPOT tells a user who plans a study                  //
// (Figure 3, Tables S9 and S10, Figures S5 and S6)                                       //
//                                                                                        //
// Seven parts, chosen by a letter on the command line:                                   //
//   A : museum insect, power against s for several starting frequencies                  //
//   B : museum insect, where to add samples (historical, modern, depth)                  //
//   C : museum insect, power against Ne for 1, 2 and 5 times the sample sizes            //
//   D : long-lived mammal, power with and without a bottleneck between dates             //
//   E : bacterial pathogen (haploid, many generations, rescaled), power against s        //
//   F : human ancient DNA, power against s for several depths of the ancient samples     //
//   G : human ancient DNA, the same number of reads per ancient date spread over         //
//       more or fewer individuals                                                        //
//                                                                                        //
// Run from the repository root, for example:                                             //
//   node paper/worked_examples/worked_examples.js A > paper/worked_examples/bee_s_p0.tsv //
////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////
// Load the simulation code //
//////////////////////////////
const S = require('../../sim.js');


/////////////////////////////////////////////////////////////
// Settings of the museum insect case study (haplodiploid) //
/////////////////////////////////////////////////////////////
const BEE = {system: 'haplodiploid',         // Inheritance system
  neByBlock: false, neBlocks: null,          // Ne is not set separately for each interval
  genTime: 1,                                // Years per generation
  Nf: 1000, Nm: 1000,                        // Breeding females and males per generation (Ne = 1,500)
  noDrift: false,                            // Drift is simulated
  p0: 0.01,                                  // Starting frequency of the favoured allele (changed in each part)
  s: 0.05,                                   // Selection coefficient (changed in each part)
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
  reps: 1000,                                // Selection simulations per run
  seed: 11};                                 // Random seed


//////////////////////////////////////////////////////////////////////
// Settings of the three other case studies                         //
// Each one = the settings above with some values changed (the same //
// values as the example scenarios of the application)              //
//////////////////////////////////////////////////////////////////////
// same values as the example scenarios of the page (PRESETS in parts/ui.js)
const MAMMAL = Object.assign({}, BEE, {system: 'diploid', genTime: 20, Nf: 2500, Nm: 2500, neByBlock: true, neBlocks: [5000, 200],
  p0: 0.1, s: 0.2, onset: 1880, years: [1880, 1950, 2020], n: [20, 20, 40], depthEarly: 4, depthLast: 15});   // Long-lived mammal: 20 years per generation, Ne 5,000 then 200
const PATHOGEN = Object.assign({}, BEE, {system: 'haploid', genTime: 0.01, Nf: 100000, Nm: 100000, p0: 0.05, s: 0.0002,
  onset: 1850, years: [1850, 1900, 1950, 2000, 2020], n: [10, 10, 15, 20, 40], depthEarly: 10, depthLast: 30});   // Bacterial pathogen: haploid, 100 generations per year (17,000 generations)
const HUMAN = Object.assign({}, BEE, {system: 'diploid', genTime: 28, Nf: 5000, Nm: 5000, p0: 0.1, s: 0.005, onset: -7000,
  years: [-7000, -5000, -3000, -1000, 1000, 2000], n: [50, 50, 50, 50, 50, 100], depthEarly: 2, depthLast: 20, err: 0.01, minDP: 1});   // Human ancient DNA: 28 years per generation, six dates over 9,000 years


////////////////////////////////////////////////////////////////////////
// Functions: one SPOT run, returning the columns of the results file //
// P: settings; extra: settings changed for this run                  //
////////////////////////////////////////////////////////////////////////
// Version used for the human case study (parts F and G): includes the quality of the genotypes
const runQ = (P, extra) => {
  const r = S.runPower(Object.assign({}, P, extra));       // One SPOT run
  return [r.trend.powerEmp.toFixed(4),                     // Power of the trend test, drift-aware threshold
    r.trend.fprNaive.toFixed(4),                           // False-positive rate of the sampling-only threshold
    r.two.powerEmp.toFixed(4),                             // Power of the first-versus-last test
    r.missRate.toFixed(4),                                 // Proportion of individuals without a genotype
    r.hetAcc.toFixed(4)];                                  // Proportion of heterozygotes called correctly
};
const HQ = ['power_trend', 'fpr_naive', 'power_first_vs_last', 'share_missing', 'het_called_correctly'];   // Names of these columns

// Version used for the other case studies (parts A to E)
const run = (P, extra) => {
  const r = S.runPower(Object.assign({}, P, extra));       // One SPOT run
  return [r.trend.powerEmp.toFixed(4),                     // Power of the trend test, drift-aware threshold
    r.trend.powerNaive.toFixed(4),                         // Power of the trend test, sampling-only threshold
    r.trend.fprNaive.toFixed(4),                           // False-positive rate of the sampling-only threshold
    r.two.powerEmp.toFixed(4)];                            // Power of the first-versus-last test
};
const H = ['power_trend', 'power_naive', 'fpr_naive', 'power_first_vs_last'];   // Names of these columns


/////////////////////////////////////
// Part chosen on the command line //
/////////////////////////////////////
const part = process.argv[2];                // Letter of the part to run
const t0 = Date.now(),                       // Starting time
  log = s => console.error(s, ((Date.now() - t0) / 1000).toFixed(0) + 's');   // Function: progress message with the time used


////////////////////////////////////////////////////////////////////
// Part A (Figure 3a): museum insect, power against the selection //
// coefficient, for four starting frequencies                     //
////////////////////////////////////////////////////////////////////
if (part === 'A') {
  console.log(['p0', 's', ...H].join('\t'));                                   // Header
  for (const p0 of [0.01, 0.05, 0.1, 0.3]) {                                   // Loop over the starting frequencies
    for (const s of [0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2])    // Loop over the selection coefficients
    console.log([p0, s, ...run(BEE, {p0, s})].join('\t'));
    log('p0 ' + p0);
  }
}


/////////////////////////////////////////////////////////////////////
// Part B (Figure 3d): museum insect, where to add effort          //
// Each option adds 60 individuals or doubles the historical depth //
/////////////////////////////////////////////////////////////////////
if (part === 'B') { // each option adds 60 individuals or doubles historical depth; 5,000 replicates to keep Monte Carlo error below 1 point
  // two regimes: Ne = 1,500 (drift strong over 90 generations) and Ne = 100,000 (drift weak, sampling error dominates);
  // in each, s is chosen so that the reference design has about 55-65% power
  const D = [['reference: 20/20/20/40', {}],
    ['+60 historical: 40/40/40/40', {n: [40, 40, 40, 40]}],
    ['+60 modern: 20/20/20/100', {n: [20, 20, 20, 100]}],
    ['historical depth 20x (instead of 10x)', {depthEarly: 20}],
    ['first and last only: 50/50 (same total)', {years: [1935, 2025], n: [50, 50]}]];   // The designs: [label, settings changed]
  const R = [[1500, 1000, [[0.01, 0.05], [0.1, 0.02], [0.3, 0.015]]],
    [100000, 66667, [[0.01, 0.035], [0.1, 0.015], [0.3, 0.012]]]];             // Two levels of drift: [Ne, breeding females (= males), list of (p0, s)]
  console.log(['Ne', 'design', 'p0', 's', 'individuals', ...H].join('\t'));    // Header
  for (const [Ne, N, cases] of R)                                              // Loop over the two levels of drift
    for (const [p0, s] of cases) {                                             // Loop over the (p0, s) pairs
      for (const [lab, o] of D) {                                              // Loop over the designs
    const P = Object.assign({}, BEE, {Nf: N, Nm: N, reps: 5000}, o);           // Settings of this design (5,000 selection simulations)
    console.log([Ne, lab, p0, s, P.n.reduce((a, b) => a + b, 0), ...run(P, {p0, s})].join('\t'));   // One results line (with the total number of individuals)
      }
      log('B ' + Ne + ' ' + p0);
    }
}


////////////////////////////////////////////////////////////////////////
// Part C (Figure 3b, c): museum insect, power and false positives of //
// the sampling-only threshold against Ne, for 1, 2 and 5 times the   //
// sample sizes                                                       //
////////////////////////////////////////////////////////////////////////
if (part === 'C') { // power and sampling-only false positives against Ne, for 1x, 2x and 5x the sample sizes, three p0;
  // the same (p0, s) as part B at Ne = 1,500 (panel d): reference power 67%, 55% and 55%
  console.log(['p0', 's', 'N_breeding_each_sex', 'Ne', 'het_retained', 'n_scale', 'n', ...H].join('\t'));   // Header
  for (const [p0, s] of [[0.1, 0.02], [0.01, 0.05], [0.3, 0.015]])             // Loop over the (p0, s) pairs
    for (const k of [1, 2, 5]) {                                               // Loop over the multipliers of the sample sizes
      for (const N of [100, 200, 300, 500, 1000, 2000, 5000, 20000, 100000]) { // Loop over the breeding numbers of each sex
    const n = BEE.n.map(x => x * k),                                           // Sample sizes multiplied by k
      P = Object.assign({}, BEE, {Nf: N, Nm: N, n, p0, s});                    // Settings of this run
    console.log([p0, s, N, Math.round(S.neFor(P)), S.hetRetained(P).toFixed(4), k, n.join('/'), ...run(P)].join('\t'));   // One results line (with Ne and the diversity left)
      }
      log('C p0 ' + p0 + ' x' + k);
    }
}


////////////////////////////////////////////////////////////
// Part D: long-lived mammal, three demographic histories //
////////////////////////////////////////////////////////////
if (part === 'D') {
  console.log(['history', 's', ...H].join('\t'));                              // Header
  for (const [lab, o] of [['bottleneck: Ne 5,000 then 200', {}],
    ['constant Ne 5,000', {neBlocks: [5000, 5000]}],
    ['severe bottleneck: Ne 5,000 then 50', {neBlocks: [5000, 50]}]]) {        // Loop over the histories: [label, settings changed]
    for (const s of [0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5])                    // Loop over the selection coefficients
      console.log([lab, s, ...run(MAMMAL, Object.assign({s}, o))].join('\t'));
    log(lab);
  }
}


///////////////////////////////////////////////////////////////
// Part E (Figure S5): bacterial pathogen, power against the //
// selection coefficient (5,000 selection simulations)       //
///////////////////////////////////////////////////////////////
if (part === 'E') {
  console.log(['s_per_generation', 's_per_year', ...H].join('\t'));            // Header
  for (const s of [0.00002, 0.00005, 0.0001, 0.0002, 0.0005]) {                // Loop over the selection coefficients per generation
    console.log([s, (s * 100).toPrecision(2), ...run(PATHOGEN, {s, reps: 5000})].join('\t'));   // One results line (s per year = 100 x s per generation)
    log('E ' + s);
  }
}


///////////////////////////////////////////////////////////////////////
// Part F (Figure S6a): human ancient DNA, depth of the five ancient //
// samples (the present-day sample stays at 20x)                     //
///////////////////////////////////////////////////////////////////////
if (part === 'F') { // human aDNA: depth of the five ancient samples (present-day sample stays at 20x); minimum depth 1 read
  console.log(['ancient_depth', 's', ...HQ].join('\t'));                       // Header
  for (const dep of [1, 2, 5, 10, 'perfect']) {                                // Loop over the depths of the ancient samples
    for (const s of [0.001, 0.002, 0.003, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.1]) {   // Loop over the selection coefficients
    const o = dep === 'perfect' ? {perfect: true, s} : {depthEarly: dep, s};   // Settings changed: perfect genotypes, or the chosen depth
    console.log([dep, s, ...runQ(HUMAN, o)].join('\t'));
    }
    log('F ' + dep);
  }
}


//////////////////////////////////////////////////////////////////////
// Part G (Figure S6b, Table S10): human ancient DNA, 100 reads per //
// ancient date spread over more or fewer individuals               //
//////////////////////////////////////////////////////////////////////
if (part === 'G') { // human aDNA: 100 reads per ancient date, spread over more or fewer individuals (5,000 replicates)
  console.log(['individuals_per_ancient_date', 'ancient_depth', 's', ...HQ].join('\t'));   // Header
  for (const s of [0.005, 0.0075])                                             // Loop over the selection coefficients
    for (const [n, dep] of [[100, 1], [50, 2], [25, 4], [20, 5], [10, 10]]) {  // Loop over (individuals per ancient date, depth): 100 reads in total
    const P = Object.assign({}, HUMAN, {n: [n, n, n, n, n, 100], depthEarly: dep, s, reps: 5000});   // Settings of this run
    console.log([n, dep, s, ...runQ(P)].join('\t'));
    log('G ' + n + ' ' + s);
  }
}
