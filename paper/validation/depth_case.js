//////////////////////////////////////////////////////////////////////////////////////
// Effect of sequencing depth in the museum insect case study (Figure 3e, Table S9) //
//                                                                                  //
// The reference case is sequenced at different depths, from 1x to 30x and with     //
// perfect genotypes. Other designs add a higher error rate per read at low depth   //
// (for example post-mortem DNA damage), unequal depth between dates, or a minimum  //
// depth of 3 reads. Everything else is the reference case.                         //
// 5,000 selection and 10,000 neutral simulations per design.                       //
//                                                                                  //
// Run from the repository root:                                                    //
//   node paper/validation/depth_case.js > paper/worked_examples/depth_case.tsv     //
//////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////
// Load the simulation code //
//////////////////////////////
const {S, with: W} = require('./base.js');   // S: simulation code; W: reference case with some settings changed


/////////////////////
// Values explored //
/////////////////////
const NE = [[1000, 1500], [66667, 100000]];   // Two levels of drift: [breeding females (= males), resulting Ne] for a haplodiploid
const SS = [0.015, 0.02, 0.03];               // Selection coefficients


////////////////////////////////////////////////////////////////////////
// List of sequencing designs                                         //
// Each design: [label, depth of the earlier dates, depth of the last //
// date, error rate per read, perfect genotypes?, minimum depth]      //
////////////////////////////////////////////////////////////////////////
const D = [];
D.push(['perfect', 10, 10, 0.005, true, 1]);                                   // Perfect genotypes (no sequencing error)
for (const d of [1, 2, 5, 10, 30])                                             // Same depth at all dates, from 1x to 30x
  D.push([`${d}x all dates`, d, d, 0.005, false, 1]);
for (const e of [0.02, 0.05])                                                  // Higher error rates per read ...
  for (const d of [1, 2])                                                      // ... at low depth (1x and 2x)
    D.push([`${d}x all dates, error ${e}`, d, d, e, false, 1]);
D.push(['2x earlier, 30x last', 2, 30, 0.005, false, 1]);                      // Low depth for the earlier dates only
D.push(['30x earlier, 2x last', 30, 2, 0.005, false, 1]);                      // Low depth for the last date only
for (const d of [1, 2])                                                        // Low depth with a minimum of 3 reads to call a genotype
  D.push([`${d}x all dates, minimum depth 3`, d, d, 0.005, false, 3]);


////////////////////////////////
// Header of the results file //
////////////////////////////////
const H = ['Ne','s','design','depth_earlier','depth_last','error','min_depth','perfect','threshold','power','power_sampling_only',
  'fpr_sampling_only','missing','het_called_correctly','haploid_called_correctly','n_selection','n_neutral','seed'];
console.log(H.join('\t'));


///////////////////////////////////////////////////////////////////
// Run SPOT for every combination of Ne, s and sequencing design //
///////////////////////////////////////////////////////////////////
let seed = 1000;                                                               // Each run gets its own random seed (1001, 1002, ...)
for (const [Nf, Ne] of NE)                                                     // Loop over the two levels of drift
  for (const s of SS)                                                          // Loop over the selection coefficients
    for (const [lab, de, dl, err, perfect, minDP] of D) {                      // Loop over the sequencing designs
  seed++;
  const r = S.runPower(W({Nf, Nm: Nf, s, depthEarly: de, depthLast: dl, err, perfect, minDP, reps: 5000, seed}));   // One SPOT run
  console.log([Ne, s, lab, de, dl, err, minDP, perfect,                        // Settings of the run
    r.trend.crit.toFixed(3),                                                   // Drift-aware threshold
    r.trend.powerEmp.toFixed(4),                                               // Power with the drift-aware threshold
    r.trend.powerNaive.toFixed(4),                                             // Power with the sampling-only threshold
    r.trend.fprNaive.toFixed(4),                                               // False-positive rate of the sampling-only threshold
    r.missRate.toFixed(4),                                                     // Proportion of individuals without a genotype
    isNaN(r.hetAcc) ? '' : r.hetAcc.toFixed(4),                                // Proportion of heterozygotes called correctly
    isNaN(r.hapAcc) ? '' : r.hapAcc.toFixed(4),                                // Proportion of haploid genotypes called correctly
    5000, r.nNeutral, seed].join('\t'));                                       // Numbers of simulations and seed
}
