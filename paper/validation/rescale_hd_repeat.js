// Repeats the two haplodiploid 60,000-generation rescaling cases with the largest power difference (rescale_results.tsv),
// 5 seeds each, exact and rescaled: the difference is run-to-run variation of the threshold, not a bias.
// Run from the repository root: node paper/validation/rescale_hd_repeat.js exact 0.0005 (or rescaled; s = 0.0005 or 0.0001)
const S=require('../../sim.js');
const base={system:'haplodiploid',genTime:1,Nf:1e4,Nm:1e4,noDrift:false,p0:0.1,h:0.5,onset:0,n:[20,20,20,40],perfect:false,depthEarly:10,depthLast:15,err:0.005,minDP:3,alpha:0.05,m:1,sMale:1,decline:1,mig:0,pMig:-1,kSib:1,maleFrac:0,reps:1000,years:[0,20000,40000,60000]};
const job=process.argv[2]; const s=+process.argv[3];
const out=[];for(let k=0;k<5;k++){const r=S.runPower(Object.assign({},base,{s,seed:500+k},job==='exact'?{noRescale:true}:{}));out.push([r.trend.powerEmp,r.trend.crit.toFixed(2)])}
console.log(job,s,JSON.stringify(out));
