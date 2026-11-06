// Engine check in a browser: for each example scenario of the page, read the exact numbers the page computes (its CSV
// export and its saved settings) and compare them with sim.js run on the same settings. Needs Playwright (Chromium).
// Run from the repository root: node validation/page_engine_check.js index.html sim.js
const {chromium} = require('playwright'); const path = require('path'); const fs = require('fs');
const S = require(path.resolve(process.argv[3] || 'sim.js'));
(async () => {
  const b = await chromium.launch(); const ctx = await b.newContext({acceptDownloads: true}); const pg = await ctx.newPage(); const errs = [];
  pg.on('pageerror', e => errs.push(String(e)));
  await pg.goto('file://' + path.resolve(process.argv[2]));
  const idle = () => pg.waitForFunction(() => { const t = document.getElementById('status').textContent; return / ms /.test(t) && !/Running/.test(t); }, null, {timeout: 120000});
  const grab = async (id) => { const [d] = await Promise.all([pg.waitForEvent('download'), pg.click('#' + id)]); return fs.readFileSync(await d.path(), 'utf8'); };
  const rows = [];
  for (const pre of ['defaults', 'mammal', 'human', 'pathogen', 'large', 'analytic']) {
    await pg.click(`[data-preset="${pre}"]`); await pg.waitForTimeout(400); await idle();
    const st = JSON.parse(await grab('save-settings')).settings;
    const csv = Object.fromEntries((await grab('res-csv')).trim().split('\n').map(l => { const i = l.indexOf(','); return [l.slice(0, i), l.slice(i + 1)]; }));
    const r = S.runPower(st);
    const g = (...re) => { for (const x of re) { const k = Object.keys(csv).find(k => x.test(k)); if (k !== undefined) return csv[k]; } return NaN; };
    const page = [g(/^power_trend_drift_aware/, /^power_trend_neutral_simulation/), g(/^power_two_point_drift_aware/, /^power_first_vs_last_(neutral_simulation|drift_aware)/),
      g(/^threshold_trend$/), g(/^threshold_two_point$/, /^threshold_first_vs_last$/), g(/^naive_power$/, /^power_trend_sampling_only/), g(/diversity_lost/)].map(Number);
    const sim = [r.trend.powerEmp, r.two.powerEmp, r.trend.crit, r.two.crit, r.trend.powerNaive, 1 - S.hetRetained(st)];
    const same = page.every((v, i) => Math.abs(v - sim[i]) <= 1e-12 * Math.max(1, Math.abs(sim[i])));
    rows.push([pre, same ? 'identical' : 'DIFFERENT', page.map(v => +v.toFixed(6)).join(' '), sim.map(v => +v.toFixed(6)).join(' ')]);
  }
  console.log('preset\tresult\tpage: power trend, power first-vs-last, threshold trend, threshold first-vs-last, sampling-only power, diversity lost\tsim.js: same');
  rows.forEach(r => console.log(r.join('\t')));
  console.log('page errors:', errs.length ? errs : 'none');
  await b.close();
})();
