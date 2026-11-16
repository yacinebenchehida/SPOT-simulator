//////////////////////////////////////////////////////////////////////////////////////
// Rescaling: repeat of the two haplodiploid cases with the largest difference (S6) //
//                                                                                  //
// In rescale_results.tsv, two haplodiploid cases of 60,000 generations showed the  //
// largest difference in power between the simulation without rescaling ("exact")   //
// and the rescaled simulation. Here each is repeated with 5 seeds, exact and       //
// rescaled. Result: the difference comes from the variation of the threshold       //
// between runs, not from a bias of the rescaling.                                  //
//                                                                                  //
// Run from the repository root (s = 0.0005 or 0.0001):                             //
//   node paper/validation/rescale_hd_repeat.js exact 0.0005                        //
//   node paper/validation/rescale_hd_repeat.js rescaled 0.0005                     //
//////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////
// Load the simulation code //
//////////////////////////////
const S=require('../../sim.js');


/////////////////////////////////
// Settings shared by all runs //
/////////////////////////////////
const base={system:'haplodiploid',   // Inheritance system
  genTime:1,                         // Years per generation
  Nf:1e4,Nm:1e4,                     // Breeding females and males per generation
  noDrift:false,                     // Drift is simulated
  p0:0.1,                            // Starting frequency of the favoured allele
  h:0.5,                             // Dominance
  onset:0,                           // Selection starts at the first date
  n:[20,20,20,40],                   // Individuals sampled at each date
  perfect:false,                     // Genotypes are called from sequencing reads
  depthEarly:10,depthLast:15,        // Mean sequencing depth (earlier dates, last date)
  err:0.005,                         // Error rate per read
  minDP:3,                           // Minimum number of reads to call a genotype
  alpha:0.05,                        // Significance level
  m:1,                               // One candidate locus
  sMale:1,decline:1,mig:0,pMig:-1,   // No sex difference in selection, constant Ne, no migration
  kSib:1,maleFrac:0,                 // Unrelated individuals, females only
  reps:1000,                         // Selection simulations per run
  years:[0,20000,40000,60000]};      // Sampling dates (one generation per year: 60,000 generations)


///////////////////////////////////
// Arguments of the command line //
///////////////////////////////////
const job=process.argv[2];           // "exact" (no rescaling) or "rescaled"
const s=+process.argv[3];            // Selection coefficient


////////////////////////////////////
// Five runs with different seeds //
////////////////////////////////////
const out=[];                        // One [power, threshold] pair per run
for(let k=0;k<5;k++){
  const r=S.runPower(Object.assign({},base,{s,seed:500+k},job==='exact'?{noRescale:true}:{}));   // One SPOT run (seeds 500 to 504), with rescaling switched off if "exact"
  out.push([r.trend.powerEmp,r.trend.crit.toFixed(2)])                                           // Keep the power and the threshold of this run
}
console.log(job,s,JSON.stringify(out));   // Print the five results on one line
