// Figure 1 of the paper: the SPOT interface for the museum insect case study (reference case, p0 = 0.1, s = 0.02).
// Opens index.html with the reference settings, draws the power curve (individuals at each earlier point) and the
// sample-size planner (selection coefficient), and saves four panels:
//   a, result tiles (screenshot); b, distribution of Z; c, power curve; d, sample-size planner
// (b-d from the page's own chart export, light colours, with titles and legends), then assembles Figure1.png.
// Needs Playwright (Chromium) and Python with Pillow. Run from the repository root: node figures/make_figure1.js
const {chromium} = require('playwright'); const path = require('path'); const fs = require('fs'); const {execFileSync} = require('child_process');
const {REF} = require('../validation/base.js');
const OUT = path.resolve('figures'); const S = 3;   // pixel density
(async () => {
  const st = Object.assign({}, REF, {errLast: REF.err});
  const hash = 's=' + Buffer.from(JSON.stringify(st)).toString('base64').replace(/\+/g, '-').replace(/\//g, '_').replace(/=+$/, '');
  const b = await chromium.launch(); const ctx = await b.newContext({viewport: {width: 1180, height: 3000}, deviceScaleFactor: S, colorScheme: 'light'});
  const p = await ctx.newPage(); const errs = []; p.on('pageerror', e => errs.push(String(e)));
  await p.goto('file://' + path.resolve('index.html') + '#' + hash);
  await p.waitForFunction(() => / ms /.test(document.getElementById('status').textContent), null, {timeout: 60000});
  await p.waitForTimeout(800);
  await p.selectOption('#sweepvar', 'nHist'); await p.click('#sweep-run');
  await p.waitForFunction(() => /settings, each/.test(document.getElementById('sweep-status').textContent), null, {timeout: 300000});
  await p.selectOption('#plan-var', 's'); await p.click('#plan-run');
  await p.waitForFunction(() => /settings, each/.test(document.getElementById('plan-status').textContent), null, {timeout: 300000});
  await (await p.$('.tiles')).screenshot({path: path.join(OUT, 'Figure1a.png')});
  for (const [id, f] of [['ch-dist', 'Figure1b'], ['ch-sweep', 'Figure1c'], ['ch-plan', 'Figure1d']]) {
    const svg = await p.evaluate(i => window.__exportChart(document.getElementById(i)).svg, id);
    fs.writeFileSync(path.join(OUT, f + '.svg'), svg);
    const q = await ctx.newPage(); await q.setContent('<html><body style="margin:0;background:#fff">' + svg + '</body></html>');
    await (await q.$('svg')).screenshot({path: path.join(OUT, f + '.png')}); await q.close();
  }
  const vals = await p.evaluate(() => ({power: document.getElementById('t-pe').textContent, status: document.getElementById('status').textContent}));
  console.log('page:', vals.power, '|', vals.status); console.log('errors:', errs.length ? errs : 'none');
  await b.close();
  execFileSync('python3', [path.join(__dirname, 'assemble_figure1.py'), OUT], {stdio: 'inherit'});
})();
