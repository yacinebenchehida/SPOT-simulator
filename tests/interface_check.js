///////////////////////////////////////////////////////////////////////////////////////
// Test: touching the interface does not change the settings used by the simulations //
//                                                                                   //
// For each example scenario, and for settings opened from a link (dates before      //
// 1935, values between slider steps, Ne per interval, the optional settings, ...),  //
// the test saves the settings, changes one setting and puts it back, and checks     //
// that the saved settings are the same. It also checks that settings opened from a  //
// link are applied as given, and that the power curve and the sample-size planner   //
// run without page errors.                                                          //
//                                                                                   //
// Needs Playwright (Chromium).                                                      //
// Run from the repository root:                                                     //
//   node tests/interface_check.js index.html                                        //
///////////////////////////////////////////////////////////////////////////////////////
const {chromium}=require('playwright'); const fs=require('fs');
const FILE=require('path').resolve(process.argv[2]||'index.html');   // The application (first argument)
(async()=>{const b=await chromium.launch();const ctx=await b.newContext({acceptDownloads:true});const p=await ctx.newPage();const errs=[];p.on('pageerror',e=>errs.push(String(e)));
const idle=()=>p.waitForFunction(()=>{const t=document.getElementById('status').textContent;return / ms /.test(t)&&!/Running/.test(t)},null,{timeout:180000});
const grab=async id=>{const [d]=await Promise.all([p.waitForEvent('download'),p.click('#'+id)]);return fs.readFileSync(await d.path(),'utf8')};
const settings=async()=>JSON.parse(await grab('save-settings')).settings;
const poke=async()=>{ // change the significance level and put it back: fires the input handler twice
  const a=await p.$eval('#alpha',e=>e.value); const other=a==='0.05'?'0.01':'0.05';
  await p.selectOption('#alpha',other); await p.waitForTimeout(150); await p.selectOption('#alpha',a); await p.waitForTimeout(300); await idle(); };
const diff=(A,B)=>Object.keys(A).filter(k=>JSON.stringify(A[k])!==JSON.stringify(B[k])).map(k=>`${k}: ${JSON.stringify(A[k])} -> ${JSON.stringify(B[k])}`);
const enc=o=>Buffer.from(unescape(encodeURIComponent(JSON.stringify(o)))).toString('base64').replace(/\+/g,'-').replace(/\//g,'_').replace(/=+$/,'');
await p.goto('file://'+FILE); await p.waitForTimeout(1200); await idle();
const report=[];
// 1. every example scenario
for(const pre of await p.$$eval('[data-preset]',a=>a.map(x=>x.dataset.preset))){
  await p.click(`[data-preset="${pre}"]`); await p.waitForTimeout(400); await idle();
  const A=await settings(); await poke(); const B=await settings(); const d=diff(A,B);
  report.push(['preset '+pre, d.length?'CHANGED '+d.join('; '):'unchanged']);
}
// 2. settings loaded from a link, including values off the slider steps and the new options
const base=await settings();
const CASES={
 'early years, late onset':{years:[1850,1900,1950,2020],n:[10,10,10,40],onset:1920},
 'BCE years':{years:[-5000,-2000,0,2000],n:[20,20,20,40],onset:-3000,genTime:25},
 'off-step values':{mig:0.001,sMale:0.33,maleFrac:0.27,depthEarly:2.3,depthLast:7.7,decline:0.37,p0:0.137,s:0.0123,Nf:1234,Nm:987},
 'Ne per interval':{neByBlock:true,neBlocks:[300,5000,800],years:[1935,1955,1975,2025],n:[20,20,20,40]},
 'new options':{calls:'pseudo',depthCV:0.5,errByDate:[0.05,0.02,0.01,0.005],neDirect:true,minDP:1},
 'haplodiploid with males and siblings':{system:'haplodiploid',maleFrac:0.4,kSib:3,sMale:0.5},
 'haploid':{system:'haploid',Nf:20000,Nm:20000},
 'no drift':{noDrift:true},
 'many loci':{m:500,alpha:0.05},
};
for(const [name,v] of Object.entries(CASES)){
  const S0=Object.assign({},base,v);
  await p.goto('about:blank'); await p.goto('file://'+FILE+'#s='+enc(S0)); await p.waitForTimeout(1200); await idle();
  const A=await settings(); const d0=diff(S0,A);
  await poke(); const B=await settings(); const d=diff(A,B);
  report.push(['link: '+name, (d0.length?'LOADED DIFFERENTLY '+d0.join('; ')+' | ':'')+(d.length?'CHANGED AFTER TOUCH '+d.join('; '):'unchanged after touch')]);
}
// 3. power curve and planner with the new options
await p.goto('about:blank'); await p.goto('file://'+FILE+'#s='+enc(Object.assign({},base,CASES['new options']))); await p.waitForTimeout(1200); await idle();
await p.click('#sweep-run'); await p.waitForTimeout(500);
await p.waitForFunction(()=>!/Running|running|…/.test(document.getElementById('sweep-run').textContent||''),null,{timeout:300000}).catch(()=>{});
await p.click('#plan-run'); await p.waitForTimeout(3000);
report.push(['curve and planner with the new options', errs.length?'ERRORS':'no page error']);
report.forEach(r=>console.log(r.join(' :: ')));
console.log('page errors:',errs.length?errs:'none'); await b.close();})();
