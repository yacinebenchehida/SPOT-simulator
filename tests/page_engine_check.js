//////////////////////////////////////////////////////////////////////////////////////
// Test: the application gives the same numbers as sim.js, in a real browser        //
//                                                                                  //
// For each example scenario of the application (and with the optional settings), //
// the test reads the exact numbers                                                 //
// the page computes (its CSV export and its saved settings) and compares them with //
// sim.js run on the same settings.                                                 //
//                                                                                  //
// Needs Playwright (Chromium).                                                     //
// Run from the repository root:                                                    //
//   node tests/page_engine_check.js index.html sim.js                              //
//////////////////////////////////////////////////////////////////////////////////////


///////////////////////////
// Load the dependencies //
///////////////////////////
const {chromium} = require('playwright');    // Browser that opens the application
const path = require('path');                // To build file names
const fs = require('fs');                    // To read the downloaded files
const S = require(path.resolve(process.argv[3] || 'sim.js'));   // The simulation code (second argument; sim.js by default)


////////////////////////////////////////////////////////////
// Open the application and compare each example scenario //
////////////////////////////////////////////////////////////
(async () => {
  // Open the page
  const b = await chromium.launch();         // Start the browser
  const ctx = await b.newContext({acceptDownloads: true});   // Allow the page to download files
  const pg = await ctx.newPage();
  const errs = [];                           // Errors of the page, if any
  pg.on('pageerror', e => errs.push(String(e)));
  await pg.goto('file://' + path.resolve(process.argv[2]));   // The application (first argument)

  // Function: wait until the page has finished its calculation
  const idle = () => pg.waitForFunction(() => {
    const t = document.getElementById('status').textContent;
    return / ms /.test(t) && !/Running/.test(t);
  }, null, {timeout: 120000});

  // Function: click a download button and return the content of the downloaded file
  const grab = async (id) => {
    const [d] = await Promise.all([pg.waitForEvent('download'), pg.click('#' + id)]);
    return fs.readFileSync(await d.path(), 'utf8');
  };

  const rows = [];                           // One results line per example scenario
  // Example scenarios, then the optional settings (version 1.1) applied to one of them
  const set = (id, v) => pg.$eval('#' + id, (e, v) => { if (e.type === 'checkbox') e.checked = v; else e.value = v; e.dispatchEvent(new Event('input', {bubbles: true})); }, v);
  const SCEN = ['defaults', 'mammal', 'human', 'pathogen', 'large', 'analytic'].map(p => [p, p, null]).concat([
    ['human + pseudo-haploid calls', 'human', async () => set('calls', 'pseudo')],
    ['human + variable depth (CV 0.5)', 'human', async () => set('dcv', '0.5')],
    ['defaults + error rate per point', 'defaults', async () => { await set('errdate', true); await pg.waitForTimeout(200); await set('er0', '0.05'); }],
    ['defaults + Ne entered directly', 'defaults', async () => { await set('nedirect', true); await pg.waitForTimeout(200); await set('ne', '40'); }]]);
  for (const [pre, preset, apply] of SCEN) {   // Loop over the scenarios
    // Load the scenario, apply the optional settings and wait for the results
    await pg.click(`[data-preset="${preset}"]`);
    await pg.waitForTimeout(400);
    await idle();
    if (apply) { await apply(); await pg.waitForTimeout(400); await idle(); }

    // Read the settings and the results of the page
    const st = JSON.parse(await grab('save-settings')).settings;   // Settings saved by the page
    const csv = Object.fromEntries((await grab('res-csv')).trim().split('\n').map(l => {
      const i = l.indexOf(',');
      return [l.slice(0, i), l.slice(i + 1)];
    }));                                     // Results of the page: name -> value

    // Run sim.js on the same settings
    const r = S.runPower(st);

    // Function: value of the first CSV line whose name matches one of the patterns (names changed between versions)
    const g = (...re) => {
      for (const x of re) {
        const k = Object.keys(csv).find(k => x.test(k));
        if (k !== undefined) return csv[k];
      }
      return NaN;
    };

    // The six values compared: from the page, then from sim.js
    const page = [g(/^power_trend_drift_aware/, /^power_trend_neutral_simulation/), g(/^power_two_point_drift_aware/, /^power_first_vs_last_(neutral_simulation|drift_aware)/),
      g(/^threshold_trend$/), g(/^threshold_two_point$/, /^threshold_first_vs_last$/), g(/^naive_power$/, /^power_trend_sampling_only/), g(/diversity_lost/)].map(Number);
    const sim = [r.trend.powerEmp, r.two.powerEmp, r.trend.crit, r.two.crit, r.trend.powerNaive, 1 - S.hetRetained(st)];
    const same = page.every((v, i) => Math.abs(v - sim[i]) <= 1e-12 * Math.max(1, Math.abs(sim[i])));   // Identical up to the rounding of the computer
    rows.push([pre, same ? 'identical' : 'DIFFERENT', page.map(v => +v.toFixed(6)).join(' '), sim.map(v => +v.toFixed(6)).join(' ')]);
  }

  // Print the comparison
  console.log('preset\tresult\tpage: power trend, power first-vs-last, threshold trend, threshold first-vs-last, sampling-only power, diversity lost\tsim.js: same');
  rows.forEach(r => console.log(r.join('\t')));
  console.log('page errors:', errs.length ? errs : 'none');
  await b.close();                           // Stop the browser
})();
