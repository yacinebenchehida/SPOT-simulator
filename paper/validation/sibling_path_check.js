/////////////////////////////////////////////////////////////////////////////////////////////
// Check of the sampling of families (siblings) (S3)                                       //
//                                                                                         //
// SPOT has two ways of sampling individuals:                                              //
//   default : each individual is drawn from the genotype frequencies of the               //
//             offspring after selection;                                                  //
//   family  : a mother and a father are drawn from the female and male offspring          //
//             after selection, and the sampled individual is their offspring              //
//             (Mendelian inheritance, no further selection).                              //
// With families of one individual (k = 1), both give unrelated individuals, so            //
// the family procedure (P.famPath in sim.js) must give the same power as the              //
// default one. This script checks it.                                                     //
//                                                                                         //
// Settings: diploid and haplodiploid; Ne = 300, 1,500 and 100,000; s = 0.02 and           //
// 0.05; design of the reference case (p0 = 0.1, 4 sampling dates, 20/20/20/40             //
// females); perfect genotypes; 20,000 selection and 40,000 neutral simulations            //
// per run; 3 seeds. Each run sets its own threshold, as SPOT does.                        //
//                                                                                         //
// se_difference: binomial standard error of the difference between two powers             //
// from 20,000 selection simulations each (it ignores the variation of the                 //
// threshold between runs, so it slightly underestimates the true standard error).         //
//                                                                                         //
// Run from the repository root (about 5 minutes):                                         //
//   node paper/validation/sibling_path_check.js > paper/validation/sibling_path_check.tsv //
/////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////
// Load the simulation code //
//////////////////////////////
const {S, with: W} = require('./base.js');   // S: simulation code; W: reference case with some settings changed


/////////////////////////////////
// Fixed values and small tool //
/////////////////////////////////
const R = 20000,                             // Selection simulations per run
  SEEDS = [1, 2, 3];                         // Random seeds
const sizes = (system, Ne) => system === 'diploid' ? {Nf: Ne / 2, Nm: Ne / 2} : {Nf: Ne / 1.5, Nm: Ne / 1.5};   // Breeding females and males that give the wanted Ne (diploid: Ne/2 each; haplodiploid: Ne/1.5 each)


///////////////////////////////////////////////////////////////
// Compare the two sampling procedures for every combination //
///////////////////////////////////////////////////////////////
console.log(['system', 'Ne', 's', 'seed', 'statistic', 'power_default', 'power_family_path', 'difference', 'se_difference', 'z'].join('\t'));   // Header
const t0 = Date.now();                                                  // Starting time, to print the time used
for (const system of ['diploid', 'haplodiploid'])                       // Loop over the inheritance systems
  for (const Ne of [300, 1500, 100000])                                 // Loop over the effective population sizes
    for (const s of [0.02, 0.05])                                       // Loop over the selection coefficients
      for (const seed of SEEDS) {                                       // Loop over the seeds
  const P = W(Object.assign({system, s, perfect: true, reps: R, seed}, sizes(system, Ne)));   // Settings of this run
  const a = S.runPower(P),                                              // Default sampling procedure
    b = S.runPower(Object.assign({}, P, {famPath: true}));              // Family sampling procedure
  for (const k of ['trend', 'two']) {                                   // For the trend test and the first-versus-last test
    const p1 = a[k].powerEmp,                                           // Power with the default procedure
      p2 = b[k].powerEmp,                                               // Power with the family procedure
      d = p2 - p1,                                                      // Difference
      se = Math.sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / R);              // Binomial standard error of the difference
    console.log([system, Ne, s, seed, k, p1.toFixed(4), p2.toFixed(4), d.toFixed(4), se.toFixed(4), se > 0 ? (d / se).toFixed(2) : ''].join('\t'));   // One results line (z = difference / standard error)
  }
  console.error(system, Ne, s, seed, ((Date.now() - t0) / 1000).toFixed(0) + 's');   // Progress message with the time used
}
