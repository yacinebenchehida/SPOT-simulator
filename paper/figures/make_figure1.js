////////////////////////////////////////////////////////////////////////////////////
// Figure 1 of the article: the SPOT interface for the museum insect case study   //
// (reference case, p0 = 0.1, s = 0.02)                                           //
//                                                                                //
// The script opens index.html with the settings of the reference case, draws the //
// power curve (individuals at each earlier point) and the sample-size planner    //
// (selection coefficient), and saves five panels:                                //
//   a : result tiles (screenshot)                                                //
//   b : distribution of Z                                                        //
//   c : power curve                                                              //
//   d : sample-size planner                                                      //
//   e : allele trajectories                                                      //
// Panels b to e come from the chart export of the application itself (light      //
// colours, with titles and legends). The script then assembles Figure1.png with  //
// assemble_figure1.py.                                                           //
//                                                                                //
// Needs Playwright (Chromium) and Python with Pillow.                            //
// Run from the repository root:                                                  //
//   node paper/figures/make_figure1.js                                           //
////////////////////////////////////////////////////////////////////////////////////


///////////////////////////
// Load the dependencies //
///////////////////////////
const {chromium} = require('playwright');    // Browser that opens the application
const path = require('path');                // To build file names
const fs = require('fs');                    // To write the SVG files
const {execFileSync} = require('child_process');   // To run the Python script that assembles the figure
const {REF} = require('../validation/base.js');    // Settings of the reference case


//////////////////
// Fixed values //
//////////////////
const OUT = path.resolve('paper/figures');   // Output folder
const S = 3;   // pixel density


////////////////////////////////////////////////////////////////////
// Open the application, run the calculations and save the panels //
////////////////////////////////////////////////////////////////////
(async () => {
  // Address of the page with the settings of the reference case (the same format as "Copy link" in the application)
  const st = Object.assign({}, REF, {errLast: REF.err});   // Settings (same error rate at the last date)
  const hash = 's=' + Buffer.from(JSON.stringify(st)).toString('base64').replace(/\+/g, '-').replace(/\//g, '_').replace(/=+$/, '');   // Settings written in the address

  // Open the page
  const b = await chromium.launch();         // Start the browser
  const ctx = await b.newContext({viewport: {width: 1180, height: 3000}, deviceScaleFactor: S, colorScheme: 'light'});   // Window size, pixel density, light colours
  const p = await ctx.newPage();
  const errs = [];                           // Errors of the page, if any
  p.on('pageerror', e => errs.push(String(e)));
  await p.goto('file://' + path.resolve('index.html') + '#' + hash);
  await p.waitForFunction(() => / ms /.test(document.getElementById('status').textContent), null, {timeout: 60000});   // Wait for the end of the main calculation
  await p.waitForTimeout(800);

  // Power curve: vary the number of individuals at each earlier point
  await p.selectOption('#sweepvar', 'nHist');
  await p.click('#sweep-run');
  await p.waitForFunction(() => /settings, each/.test(document.getElementById('sweep-status').textContent), null, {timeout: 300000});   // Wait for the curve

  // Sample-size planner: columns = selection coefficient
  await p.selectOption('#plan-var', 's');
  await p.click('#plan-run');
  await p.waitForFunction(() => /settings, each/.test(document.getElementById('plan-status').textContent), null, {timeout: 300000});   // Wait for the planner

  // Panel a: screenshot of the result tiles
  await (await p.$('.tiles')).screenshot({path: path.join(OUT, 'Figure1a.png')});

  // Panels b to e: chart export of the application, saved as SVG and as PNG
  for (const [id, f] of [['ch-dist', 'Figure1b'], ['ch-sweep', 'Figure1c'], ['ch-plan', 'Figure1d'], ['ch-traj', 'Figure1e']]) {
    const svg = await p.evaluate(i => window.__exportChart(document.getElementById(i)).svg, id);   // SVG of the chart, made by the application
    fs.writeFileSync(path.join(OUT, f + '.svg'), svg);
    const q = await ctx.newPage();           // New page that shows only this SVG ...
    await q.setContent('<html><body style="margin:0;background:#fff">' + svg + '</body></html>');
    await (await q.$('svg')).screenshot({path: path.join(OUT, f + '.png')});   // ... to save a picture of it
    await q.close();
  }

  // Print the power shown by the page and any error, then close the browser
  const vals = await p.evaluate(() => ({power: document.getElementById('t-pe').textContent, status: document.getElementById('status').textContent}));
  console.log('page:', vals.power, '|', vals.status);
  console.log('errors:', errs.length ? errs : 'none');
  await b.close();

  // Assemble the five panels into Figure1.png
  execFileSync('python3', [path.join(__dirname, 'assemble_figure1.py'), OUT], {stdio: 'inherit'});
})();
