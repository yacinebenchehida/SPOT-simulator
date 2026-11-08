// Draws the illustration at the top of the README (assets/spot_illustration.gif) from the animation of the application.
// Run from the root of the repository: node assets/make_illustration.js   (needs Playwright and ffmpeg)
// The drawing code is taken from index.html; only the motion speeds are rounded so that the animation loops without a jump.
const {chromium}=require('playwright'), fs=require('fs'), path=require('path'), {execSync}=require('child_process');
const W=900, H=420, DPR=2, FPS=20, SECONDS=6;                 // size in CSS pixels, pixel density, frames per second, loop length
const OUT=path.resolve('assets'), TMP=fs.mkdtempSync('/tmp/spot_gif_');
const html=fs.readFileSync('index.html','utf8');
let code=html.slice(html.indexOf("const cv=document.getElementById('trace')"), html.indexOf("start(); let tm"));
code=code.replace("document.getElementById('trace')","document.getElementById('trace_readme')")
  .replace('Math.sin(t*2.2)*3','Math.sin(t*2*Math.PI/3)*3')            // bee bobbing: 2 cycles in 6 s
  .replace('Math.sin(t*38)','Math.sin(t*2*Math.PI*36/6)')            // bee wings: 36 beats in 6 s
  .replace('Math.sin(t*44+1)','Math.sin(t*2*Math.PI*42/6+1)')        // hoverfly wings: 42 beats in 6 s
  .replace('Math.sin(t*1.7+1)*3','Math.sin(t*2*Math.PI/3+1)*3')      // hoverfly bobbing: 2 cycles in 6 s
  +'\nwindow.__frame=frame; window.__t0=t0; size();';
code='(function(){'+code+'})();';
(async()=>{
  const b=await chromium.launch(); const p=await b.newPage({viewport:{width:W+40,height:H+40},deviceScaleFactor:DPR,colorScheme:'light'});
  await p.addInitScript(()=>{window.requestAnimationFrame=()=>0;});
  await p.goto('file://'+path.resolve('index.html')); await p.waitForTimeout(500);
  await p.evaluate(([W,H])=>{document.body.innerHTML=`<div id="box" class="hero-art" style="position:fixed;left:0;top:0;width:${W}px;height:${H}px;min-height:0;border-radius:0"><canvas id="trace_readme"></canvas></div>`;},[W,H]);
  await p.evaluate(code);
  const box=await p.$('#box'); const n=FPS*SECONDS;
  for(let i=0;i<n;i++){
    await p.evaluate(t=>window.__frame(window.__t0+t*1000), i/FPS);
    await box.screenshot({path:path.join(TMP,String(i).padStart(4,'0')+'.png')});
  }
  await b.close();
  execSync(`ffmpeg -y -loglevel error -framerate ${FPS} -i ${TMP}/%04d.png -vf "scale=${W*DPR}:-1:flags=lanczos,split[a][b];[a]palettegen=stats_mode=diff[p];[b][p]paletteuse=dither=none:diff_mode=rectangle" -loop 0 ${OUT}/spot_illustration.gif`);
  fs.copyFileSync(path.join(TMP,'0000.png'), path.join(OUT,'spot_illustration.png'));    // still image (first frame)
  console.log('written', OUT+'/spot_illustration.gif', (fs.statSync(OUT+'/spot_illustration.gif').size/1e6).toFixed(2)+' MB');
})();
