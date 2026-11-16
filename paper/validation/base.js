////////////////////////////////////////////////////////////////////////////////////
// Reference case of the article (museum insect case study, Table 3)              //
//                                                                                //
// This file is loaded by the other scripts. It gives them:                       //
//   S    : the simulation code of SPOT (sim.js)                                  //
//   REF  : the settings of the reference case                                    //
//   with : a function that returns the reference case with some settings changed //
////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////
// Load the simulation code of SPOT //
//////////////////////////////////////
const S = require('../../sim.js');   // Same code as the one used by the application


////////////////////////////////////
// Settings of the reference case //
////////////////////////////////////
const REF = {system:'haplodiploid',  // Inheritance system
  neByBlock:false,                   // Ne is not set separately for each interval
  neBlocks:null,                     // (no Ne per interval)
  genTime:1,                         // Years per generation
  Nf:1000,                           // Breeding females per generation
  Nm:1000,                           // Breeding males per generation (Ne = 1,500 for a haplodiploid)
  noDrift:false,                     // Drift is simulated
  p0:0.1,                            // Starting frequency of the favoured allele
  s:0.02,                            // Selection coefficient
  h:0.5,                             // Dominance
  onset:1935,                        // Year selection starts
  years:[1935,1955,1975,2025],       // Sampling years
  n:[20,20,20,40],                   // Individuals sampled at each date
  perfect:false,                     // Genotypes are called from sequencing reads (not perfect)
  depthEarly:10,                     // Mean sequencing depth of the earlier dates
  depthLast:15,                      // Mean sequencing depth of the last date
  err:0.005,                         // Error rate per read
  minDP:3,                           // Minimum number of reads to call a genotype
  alpha:0.05,                        // Significance level
  m:1,                               // Number of candidate loci tested
  sMale:1,                           // Selection in males relative to females
  decline:1,                         // Ne at the last date relative to the first (1 = constant)
  mig:0,                             // Migration rate per generation
  pMig:-1,                           // Allele frequency of the migrants (-1 = same as the starting frequency)
  kSib:1,                            // Sampled individuals per family (1 = unrelated individuals)
  maleFrac:0,                        // Proportion of males in the samples
  reps:1000,                         // Number of selection simulations
  seed:1};                           // Random seed


////////////////////////////////////
// What the other scripts receive //
////////////////////////////////////
module.exports = {S,                    // The simulation code
  REF,                                  // The reference case
  with: o => Object.assign({}, REF, o)  // The reference case with the settings in "o" changed
};
