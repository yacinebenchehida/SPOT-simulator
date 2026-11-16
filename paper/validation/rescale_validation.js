//////////////////////////////////////////////////////////////////////////////////////////
// Rescaling: simulations without rescaling ("exact") against rescaled ones (S6)        //
//                                                                                      //
// For many designs, the same power calculation is run twice: once simulating every     //
// generation (exact), once with the rescaling of SPOT. The script compares power,      //
// thresholds, false-positive rates, the distribution of the allele frequency at        //
// the last date, and the run times.                                                    //
//                                                                                      //
// Run from the repository root:                                                        //
//   node paper/validation/rescale_validation.js > paper/validation/rescale_results.tsv //
//                                                                                      //
// Run times only (with nothing else running on the computer):                          //
//   node paper/validation/rescale_validation.js times [G]                              //
//   (G = 5000, 20000 or 60000; without G, all; with G, no header unless G = 5000)      //
//   -> paper/validation/rescale_times.tsv                                              //
//   It reruns the same two calculations per test case (exact, then rescaled) and       //
//   records only their run times.                                                      //
//////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////
// Load the simulation code //
//////////////////////////////
const S=require('../../sim.js');


////////////////////////////////////
// Settings shared by all designs //
////////////////////////////////////
const base={system:'haploid',        // Inheritance system (changed below)
  genTime:1,                         // Years per generation
  Nf:1e4,Nm:1e4,                     // Breeding females and males (changed below)
  noDrift:false,                     // Drift is simulated
  p0:0.1,                            // Starting frequency of the favoured allele (changed below)
  s:0.0005,                          // Selection coefficient per generation (changed below)
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
  seed:1};                           // Random seed
const designFor=G=>({years:[0,Math.round(G/3),Math.round(2*G/3),G]});   // Four sampling dates equally spaced over G generations


///////////////////////////////////////////////////////////////////////
// Function: allele frequency of the population at the last sampling //
// date, in "reps" simulations (no individuals sampled)              //
// Used to compare its distribution with and without rescaling       //
///////////////////////////////////////////////////////////////////////
// distribution of the true population frequency at the last sampling point (no sampling), exact vs rescaled
function finalFreq(P,s,reps,seed){
  const Q=S.prep(Object.assign({},P,{n:P.n.map(()=>0)})),   // Same settings, with no individual sampled
    R=S.makeRng(seed);                                     // Random number generator
  const out=new Float64Array(reps);                        // Final frequency of each simulation
  for(let r=0;r<reps;r++){
    const o=S.simulateRep(Q,s,R,true);                     // One simulation, keeping the trajectory of the population frequency
    out[r]=o.traj[o.traj.length-1];                        // Frequency at the last date
  }
  return out;
}


////////////////////////////////////////////////////////////////////////
// Function: mean, standard deviation, proportion lost and proportion //
// fixed of a list of allele frequencies                              //
////////////////////////////////////////////////////////////////////////
const stats=a=>{
  let m=0;
  for(const x of a)m+=x;
  m/=a.length;                                             // Mean
  let v=0;
  for(const x of a)v+=(x-m)*(x-m);
  v/=a.length-1;                                           // Variance
  let lost=0,fix=0;
  for(const x of a){
    if(x<=0)lost++;                                        // Allele lost
    if(x>=1)fix++                                          // Allele fixed
  }
  return {mean:m,sd:Math.sqrt(v),lost:lost/a.length,fixed:fix/a.length}
};


///////////////////////////////////
// Arguments of the command line //
///////////////////////////////////
const TIMES=process.argv[2]==='times',                     // true: record run times only
  ONLY=process.argv[3]?+process.argv[3]:null;              // Number of generations to run (null: all)


//////////////////////////////////////
// Headers of the two output tables //
//////////////////////////////////////
if(TIMES&&(!ONLY||ONLY===5000)) console.log(['G','system','N','s','p0','lambda','ms_exact','ms_rescaled','speedup'].join('\t'));   // Header of the run-time table
const rows=[['G','system','N','s','p0','lambda','scaled_s','scaled_N','power_exact','power_rescaled','two_exact','two_rescaled','crit_exact','crit_rescaled','naiveFPR_exact','naiveFPR_rescaled','meanFinal_exact','meanFinal_rescaled','sdFinal_exact','sdFinal_rescaled','neutralSD_exact','neutralSD_rescaled','neutralSD_theory','ms_exact','ms_rescaled']];   // Header of the results table


///////////////////////////////////////////////////////////////////////
// Loop over the designs: number of generations, inheritance system, //
// population size, selection coefficient and starting frequency     //
///////////////////////////////////////////////////////////////////////
const t0=Date.now();                                       // Starting time, to print the time used
for(const G of [5000,20000,60000])
  for(const system of ['haploid','diploid','haplodiploid'])
    for(const N of [1e3,1e4,1e5])
      for(const s of [0.0001,0.0005,0.002])
        for(const p0 of [0.01,0.1]){
  if(G===60000&&(N<1e4||s>0.0005||p0<0.1)) continue; // long runs: a subset, to keep the exact simulations affordable
  const P=Object.assign({},base,designFor(G),{system,Nf:N,Nm:N,s,p0});   // Settings of this design
  const sc=S.rescale(P);                                   // Rescaling chosen by SPOT (scaling factor and rescaled settings)
  if(sc.lambda===1) continue;                              // Skip the designs that SPOT does not rescale

  // Run times only: time the exact and the rescaled calculation, then go to the next design
  if(TIMES){
    if(ONLY&&G!==ONLY) continue;
    let t=Date.now();
    S.runPower(Object.assign({},P,{noRescale:true}));      // Exact calculation
    const te=Date.now()-t;                                 // Its run time (milliseconds)
    t=Date.now();
    S.runPower(Object.assign({},P,{seed:2}));              // Rescaled calculation
    const tr=Date.now()-t;                                 // Its run time (milliseconds)
    console.log([G,system,N,s,p0,sc.lambda,te,tr,(te/tr).toFixed(2)].join('\t'));
    console.error(G,system,N,s,p0,te,tr);
    continue;
  }

  // Power, thresholds and false-positive rates, exact and rescaled
  let t=Date.now();
  const e=S.runPower(Object.assign({},P,{noRescale:true}));   // Exact calculation (seed 1)
  const te=Date.now()-t;                                   // Its run time
  t=Date.now();
  const r=S.runPower(Object.assign({},P,{seed:2}));        // Rescaled calculation (seed 2)
  const tr=Date.now()-t;                                   // Its run time

  // Allele frequency at the last date: with selection (fe, fr) and without (ne, nr), exact and rescaled
  const fe=stats(finalFreq(Object.assign({},P,{noRescale:true}),s,1000,3)),
    fr=stats(finalFreq(sc.P,sc.P.s,1000,4));
  const ne=stats(finalFreq(Object.assign({},P,{noRescale:true}),0,1000,5)),
    nr=stats(finalFreq(sc.P,0,1000,6));

  // Standard deviation of the neutral frequency expected from theory
  const Ne=S.neFor(P),                                     // Effective population size
    loss=system==='haploid'?1-Math.pow(1-1/Ne,G):1-Math.pow(1-1/(2*Ne),G),   // Expected loss of diversity over G generations
    sdTheory=Math.sqrt(p0*(1-p0)*loss);                    // Expected standard deviation of the allele frequency

  // One results line per design
  rows.push([G,system,N,s,p0,sc.lambda,+sc.P.s.toPrecision(3),+sc.P.Nf.toPrecision(3),e.trend.powerEmp,r.trend.powerEmp,e.two.powerEmp,r.two.powerEmp,+e.trend.crit.toFixed(3),+r.trend.crit.toFixed(3),e.trend.fprNaive,r.trend.fprNaive,
    +fe.mean.toFixed(4),+fr.mean.toFixed(4),+fe.sd.toFixed(4),+fr.sd.toFixed(4),+ne.sd.toFixed(4),+nr.sd.toFixed(4),+sdTheory.toFixed(4),te,tr]);
  console.error(rows.length-1, G,system,N,s,p0,'lambda',sc.lambda,'power',e.trend.powerEmp,r.trend.powerEmp,((Date.now()-t0)/1000).toFixed(0)+'s');   // Progress message
}


/////////////////////////////
// Write the results table //
/////////////////////////////
if(!TIMES) console.log(rows.map(r=>r.join('\t')).join('\n'));
