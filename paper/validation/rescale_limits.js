//////////////////////////////////////////////////////////////////////////////////////
// Rescaling: stress test with forced scaling factors (S6, limits of the rescaling) //
//                                                                                  //
// Larger and larger scaling factors (lambda) are forced, and the power is compared //
// with the simulation without rescaling. This shows where the rescaling stops      //
// being accurate.                                                                  //
//                                                                                  //
// Output: one row per design and forced lambda, with the rescaled s, the smallest  //
// rescaled breeding number, the power of the trend test without rescaling          //
// (lambda = 1, seed 1) and with the forced lambda (seed 2), their difference and   //
// its binomial standard error (2,000 selection simulations each; it ignores the    //
// variation of the threshold, so it is slightly too small).                        //
//                                                                                  //
// Run from the repository root (about 1 minute):                                   //
//   node paper/validation/rescale_limits.js > paper/validation/rescale_limits.tsv  //
//////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////
// Load the simulation code //
//////////////////////////////
const S=require('../../sim.js');


//////////////////////////////////
// Settings of the first design //
//////////////////////////////////
const base={system:'haploid',        // Inheritance system
  genTime:1,                         // Years per generation
  Nf:1e4,Nm:1e4,                     // Breeding individuals (haploid: Nf is the population size)
  noDrift:false,                     // Drift is simulated
  p0:0.1,                            // Starting frequency of the favoured allele
  s:0.0005,                          // Selection coefficient per generation
  h:0.5,                             // Dominance
  onset:0,                           // Selection starts at the first date
  years:[0,6667,13333,20000],        // Sampling dates (20,000 generations)
  n:[20,20,20,40],                   // Individuals sampled at each date
  perfect:false,                     // Genotypes are called from sequencing reads
  depthEarly:10,depthLast:15,        // Mean sequencing depth (earlier dates, last date)
  err:0.005,                         // Error rate per read
  minDP:3,                           // Minimum number of reads to call a genotype
  alpha:0.05,                        // Significance level
  m:1,                               // One candidate locus
  sMale:1,decline:1,mig:0,pMig:-1,   // No sex difference in selection, constant Ne, no migration
  kSib:1,maleFrac:0,                 // Unrelated individuals, females only
  reps:2000,                         // Selection simulations per run
  seed:1};                           // Random seed


//////////////////////////////////////////////////////////////
// Function: the settings P rescaled by a forced factor lam //
// One simulated generation stands for lam real generations //
//////////////////////////////////////////////////////////////
function force(P,lam){
  const Q=Object.assign({},P,{noRescale:true});   // Copy of the settings; the automatic rescaling of SPOT is switched off
  Q.genTime=P.genTime*lam;                        // Each simulated generation lasts lam times longer
  Q.s=Math.pow(1+P.s,lam)-1;                      // Selection coefficient over lam generations
  Q.Nf=P.Nf/lam;                                  // Breeding females divided by lam (same total drift)
  Q.Nm=P.Nm/lam;                                  // Breeding males divided by lam
  return Q;
}


/////////////////////////////////////////////////////////////////
// The four designs: [label, settings that differ from "base"] //
/////////////////////////////////////////////////////////////////
const cases=[['haploid N=1e4 s=5e-4 p0=0.1',{}],
  ['haploid N=1e3 s=2e-3 p0=0.1',{Nf:1e3,Nm:1e3,s:0.002}],
  ['diploid N=1e5 s=5e-4 p0=0.01',{system:'diploid',Nf:1e5,Nm:1e5,p0:0.01}],
  ['diploid N=1e3 s=5e-4 p0=0.1',{system:'diploid',Nf:1e3,Nm:1e3}]];


//////////////////////////////////////////////////
// Compare the power with and without rescaling //
//////////////////////////////////////////////////
console.log(['design','lambda','s_rescaled','N_rescaled','power_no_rescaling','power_rescaled','difference','se_difference'].join('\t'));   // Header
for(const [lab,o] of cases){                                               // Loop over the four designs
  const P=Object.assign({},base,o);                                        // Settings of this design
  const ex=S.runPower(Object.assign({},P,{noRescale:true})).trend.powerEmp;   // Power without rescaling (reference value)
  for(const lam of [2,5,10,20,50,100,200]){                                // Loop over the forced scaling factors
    const Q=force(P,lam);                                                  // Rescaled settings
    if(Q.Nf<2) break;                                                      // Stop when the rescaled population is smaller than 2
    const r=S.runPower(Object.assign({},Q,{seed:2})).trend.powerEmp;       // Power with the forced scaling factor
    const R=P.reps,                                                        // Number of selection simulations
      se=Math.sqrt((ex*(1-ex)+r*(1-r))/R);                                 // Binomial standard error of the difference between the two powers
    console.log([lab,lam,Q.s.toPrecision(2),Math.round(Q.Nf),(ex*100).toFixed(1),(r*100).toFixed(1),((r-ex)*100).toFixed(1),(se*100).toFixed(1)].join('\t'));   // One results line (powers in %)
  }
  console.error(lab,'done');                                               // Progress message on the screen
}
