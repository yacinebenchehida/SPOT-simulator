//////////////////////////////////////////////////////////////////////////////////////
// Numbers of the explained guide (miscellaneous/simulations_explained.md)          //
//                                                                                  //
// This script reproduces the numbers of the guide that come from SPOT. Each block  //
// below changes one setting and prints the power (and other values) on the screen. //
// It uses 2,000 selection simulations and 4,000 neutral simulations (seed 1) and   //
// takes about a minute.                                                            //
//                                                                                  //
// Run from the repository root:                                                    //
//   node miscellaneous/explained_examples.js                                       //
//////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////
// Load the simulation code //
//////////////////////////////
const S=require('../sim.js');


//////////////////////////////////////////////////////////////////////
// Settings the application opens with for the "Bumblebees" example //
//////////////////////////////////////////////////////////////////////
const D={system:'haplodiploid',      // Inheritance system
  genTime:1,                         // Years per generation
  Nf:1000,Nm:1000,                   // Breeding females and males per generation
  noDrift:false,                     // Drift is simulated
  p0:0.01,                           // Starting frequency of the favoured allele
  s:0.05,                            // Selection coefficient
  h:0.5,                             // Dominance
  onset:1935,                        // Year selection starts
  years:[1935,1955,1975,2025],       // Sampling years
  n:[20,20,20,40],                   // Individuals sampled at each date
  perfect:false,                     // Genotypes are called from sequencing reads
  depthEarly:10,depthLast:15,        // Mean sequencing depth (earlier dates, last date)
  err:0.005,                         // Error rate per read
  minDP:3,                           // Minimum number of reads to call a genotype
  alpha:0.05,                        // Significance level
  m:1,                               // One candidate locus
  sMale:1,decline:1,mig:0,pMig:-1,   // No sex difference in selection, constant Ne, no migration
  kSib:1,maleFrac:0,                 // Unrelated individuals, females only
  reps:2000,                         // Selection simulations
  seed:1};                           // Random seed


/////////////////
// Small tools //
/////////////////
const R=(o)=>{                       // Function: one SPOT run with the settings in "o" changed
  const r=S.runPower({...D,...o});
  return r
};
const f=x=>(x*100).toFixed(1)+'%';   // Function: write a proportion as a percentage with one decimal


/////////////////////////////////////////////////////////////////////////
// The default example, haplodiploid then diploid                      //
// Printed: power of the trend test and of the first-versus-last test, //
// and the allele frequency at the last date without drift ("det end") //
/////////////////////////////////////////////////////////////////////////
let r=R({});
console.log('DEFAULT hap', f(r.trend.powerEmp), f(r.two.powerEmp), 'fpr',f(r.trend.fprNaive),'het',f(r.hetAcc), 'det end', S.deterministic(D).slice(-1)[0].toFixed(3), 'Ne', S.neFor(D), 'Hlost', f(1-S.hetRetained(D)));
r=R({system:'diploid'});
console.log('DEFAULT dip', f(r.trend.powerEmp), f(r.two.powerEmp), 'det end', S.deterministic({...D,system:'diploid'}).slice(-1)[0].toFixed(3));


//////////////////////////////////////////////////////////////////
// With and without drift (diploid, p0 = 0.3, s = 0.01, perfect //
// genotypes)                                                   //
//////////////////////////////////////////////////////////////////
console.log('-- no drift vs drift (diploid, p0 .3, s .01, n 20/20/20/40, perfect)');
for(const nd of [true,false]){       // Drift switched off, then on (Ne = 10,000)
  r=R({system:'diploid',p0:0.3,s:0.01,perfect:true,noDrift:nd,Nf:5000,Nm:5000});
  console.log(nd?'nodrift':'drift Ne10k', f(r.trend.powerEmp), f(r.two.powerEmp),'crit2',r.two.crit.toFixed(2),'fpr',f(r.trend.fprNaive));
}


/////////////////////////////////////////////////////////////////////
// Very large samples (500 individuals per date), with and without //
// drift                                                           //
/////////////////////////////////////////////////////////////////////
console.log('-- census-like huge n with drift (diploid p0 .3 s .01 Ne 10k)');
r=R({system:'diploid',p0:0.3,s:0.01,perfect:true,Nf:5000,Nm:5000,n:[500,500,500,500]});
console.log('n500', f(r.trend.powerEmp), f(r.two.powerEmp));
r=R({system:'diploid',p0:0.3,s:0.01,perfect:true,noDrift:true,n:[500,500,500,500]});
console.log('n500 nodrift', f(r.trend.powerEmp));


////////////////////////////////////////////////////////////////////
// False-positive rate of the sampling-only threshold against the //
// breeding numbers (haplodiploid, p0 = 0.3, s = 0.02)            //
////////////////////////////////////////////////////////////////////
console.log('-- naive FPR vs breeding numbers (hap, p0 .3, s .02)');
for(const N of [100,300,1000,3000,30000]){   // Breeding females (= males)
  r=R({p0:0.3,s:0.02,Nf:N,Nm:N});
  console.log(N, 'Ne',S.neFor({system:'haplodiploid',Nf:N,Nm:N}), 'pow',f(r.trend.powerEmp),'naivePow',f(r.trend.powerNaive),'fpr',f(r.trend.fprNaive));
}


////////////////////////////////////////////////
// Sample sizes multiplied by 0.5, 1, 2 and 4 //
////////////////////////////////////////////////
console.log('-- sample size factor (defaults hap)');
for(const k of [0.5,1,2,4]){         // Default example
  r=R({n:D.n.map(v=>Math.round(v*k))});
  console.log(k, D.n.map(v=>Math.round(v*k)).join('/'), f(r.trend.powerEmp));
}
console.log('-- sample size factor (diploid p0 .3 s .01 Ne 10k)');
for(const k of [0.5,1,2,4]){         // Diploid, weak drift
  r=R({system:'diploid',p0:0.3,s:0.01,Nf:5000,Nm:5000,n:D.n.map(v=>Math.round(v*k))});
  console.log(k, f(r.trend.powerEmp));
}


/////////////////////////////////////////////////////
// Dominance (0, 0.5, 1), haplodiploid and diploid //
/////////////////////////////////////////////////////
console.log('-- dominance hap vs dip');
for(const h of [0,0.5,1])
  for(const sys of ['haplodiploid','diploid']){
    r=R({h,system:sys});
    console.log(h,sys,f(r.trend.powerEmp),'det',S.deterministic({...D,h,system:sys}).slice(-1)[0].toFixed(4));
  }


///////////////////////////////////////////////////////////////////////
// Sequencing depth, for samples of females and for samples of males //
///////////////////////////////////////////////////////////////////////
console.log('-- depth');
for(const d of [2,3,6,10,20]){       // Mean depth at all dates
  r=R({depthEarly:d,depthLast:d});   // Females only
  console.log(d,'pow',f(r.trend.powerEmp),'het',f(r.hetAcc),'miss',f(r.missRate));
  r=R({depthEarly:d,depthLast:d,maleFrac:1});   // Males only
  console.log('  males',f(r.trend.powerEmp),'hap',f(r.hapAcc),'miss',f(r.missRate));
}


///////////////////////////////////
// One setting changed at a time //
///////////////////////////////////
console.log('-- ksib');               // Sampled individuals per family
for(const k of [1,2,5,10]){
  r=R({kSib:k});
  console.log(k,f(r.trend.powerEmp),'fpr',f(r.trend.fprNaive));
}
console.log('-- decline');            // Ne at the last date relative to the first
for(const d of [1,0.3,0.1]){
  r=R({decline:d});
  console.log(d,f(r.trend.powerEmp),'Hlost',f(1-S.hetRetained({...D,decline:d})));
}
console.log('-- migration');          // Migration rate per generation
for(const m of [0,0.01,0.05]){
  r=R({mig:m});
  console.log(m,f(r.trend.powerEmp));
}
console.log('-- onset');              // Year selection starts
for(const y of [1935,1975,2000]){
  r=R({onset:y});
  console.log(y,f(r.trend.powerEmp),f(r.two.powerEmp));
}
console.log('-- m');                  // Number of candidate loci
for(const m of [1,10,100]){
  r=R({m});
  console.log(m,f(r.trend.powerEmp),r.trend.method,r.trend.crit.toFixed(2));
}
console.log('-- sMale');              // Selection in males relative to females (0: none; 1: the same)
for(const sm of [0,1]){
  r=R({sMale:sm});
  console.log(sm,f(r.trend.powerEmp));
}
console.log('-- males sampling');     // Females, the same number of males, twice as many males (same number of gene copies)
for(const [lab,o] of [['fem',{}],['males same ind',{maleFrac:1}],['males same alleles',{maleFrac:1,n:[40,40,40,80]}]]){
  r=R(o);
  console.log(lab,f(r.trend.powerEmp));
}
console.log('-- 4pt vs 2pt late onset');   // Selection starting late: trend test against first-versus-last test
r=R({onset:1980});
console.log(f(r.trend.powerEmp),f(r.two.powerEmp));
