///////////////////////////////////////////////////////////////////////////////////
// Draw the illustration at the top of the README (assets/spot_illustration.gif) //
//                                                                               //
// The illustration is the animation shown at the top of the application. The    //
// drawing code is taken from index.html; only the speeds of the movements are   //
// rounded, so that the animation loops without a jump after 6 seconds.          //
//                                                                               //
// Needs Playwright (Chromium) and ffmpeg.                                       //
// Run from the repository root:                                                 //
//   node assets/make_illustration.js                                            //
///////////////////////////////////////////////////////////////////////////////////


///////////////////////////
// Load the dependencies //
///////////////////////////
const {chromium}=require('playwright'),      // Browser that draws the animation
  fs=require('fs'),                          // To read and copy files
  path=require('path'),                      // To build file names
  {execSync}=require('child_process');       // To run ffmpeg


//////////////////////////////////////
// Size and length of the animation //
//////////////////////////////////////
const W=900, H=520, DPR=2, FPS=20, SECONDS=6;                 // size in CSS pixels, pixel density, frames per second, loop length
const OUT=path.resolve('assets'),            // Output folder
  TMP=fs.mkdtempSync('/tmp/spot_gif_');      // Temporary folder for the frames


/////////////////////////////////////////////////////////////////
// Take the drawing code of the animation from the application //
/////////////////////////////////////////////////////////////////
const html=fs.readFileSync('index.html','utf8');   // The application
let code=html.slice(html.indexOf("const cv=document.getElementById('trace')"), html.indexOf("start(); let tm"));   // The part of the page that draws the animation
code=code.replace("document.getElementById('trace')","document.getElementById('trace_readme')")   // Draw in the canvas created below
  .replace('Math.sin(t*2.2)*3','Math.sin(t*2*Math.PI/3)*3')            // bee bobbing: 2 cycles in 6 s
  .replace('Math.sin(t*38)','Math.sin(t*2*Math.PI*36/6)')            // bee wings: 36 beats in 6 s
  .replace('Math.sin(t*44+1)','Math.sin(t*2*Math.PI*42/6+1)')        // hoverfly wings: 42 beats in 6 s
  .replace('Math.sin(t*1.7+1)*3','Math.sin(t*2*Math.PI/3+1)*3')      // hoverfly bobbing: 2 cycles in 6 s
  +'\nwindow.__frame=frame; window.__t0=t0; size();';                // Give this script access to the function that draws one frame
code='(function(){'+code+'})();';            // Wrap the code so that it can be run in the page


///////////////////////////////////////////////////////////////
// Draw the frames one by one, then assemble them into a GIF //
///////////////////////////////////////////////////////////////
(async()=>{
  // Open the application, with its own animation stopped
  const b=await chromium.launch();           // Start the browser
  const p=await b.newPage({viewport:{width:W+40,height:H+40},deviceScaleFactor:DPR,colorScheme:'light'});
  await p.addInitScript(()=>{window.requestAnimationFrame=()=>0;});   // The page must not animate by itself: this script chooses the time of each frame
  await p.goto('file://'+path.resolve('index.html'));
  await p.waitForTimeout(500);

  // Replace the page by an empty box of the wanted size, and run the drawing code in it
  await p.evaluate(([W,H])=>{document.body.innerHTML=`<div id="box" class="hero-art" style="position:fixed;left:0;top:0;width:${W}px;height:${H}px;min-height:0;border-radius:0"><canvas id="trace_readme"></canvas></div>`;},[W,H]);
  await p.evaluate(code);

  // Draw and save each frame
  const box=await p.$('#box');
  const n=FPS*SECONDS;                       // Number of frames (120)
  for(let i=0;i<n;i++){
    await p.evaluate(t=>window.__frame(window.__t0+t*1000), i/FPS);                    // Draw the animation at time i/FPS seconds
    await box.screenshot({path:path.join(TMP,String(i).padStart(4,'0')+'.png')});     // Save the frame as 0000.png, 0001.png, ...
  }
  await b.close();                           // Stop the browser

  // Assemble the frames into a GIF that loops, and keep the first frame as a still image
  execSync(`ffmpeg -y -loglevel error -framerate ${FPS} -i ${TMP}/%04d.png -vf "scale=${W*DPR}:-1:flags=lanczos,split[a][b];[a]palettegen=stats_mode=diff[p];[b][p]paletteuse=dither=none:diff_mode=rectangle" -loop 0 ${OUT}/spot_illustration.gif`);
  fs.copyFileSync(path.join(TMP,'0000.png'), path.join(OUT,'spot_illustration.png'));    // still image (first frame)
  console.log('written', OUT+'/spot_illustration.gif', (fs.statSync(OUT+'/spot_illustration.gif').size/1e6).toFixed(2)+' MB');
})();
