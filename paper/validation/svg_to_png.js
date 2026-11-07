// SVG -> PNG at 300 dpi (an SVG pixel is 1/96 inch, so the scale factor is 300/96), with the Chromium of Playwright.
// Run from the repository root: node paper/validation/svg_to_png.js file1.svg [file2.svg ...]  ->  file1.png ...
// Optional: --out DIR writes the PNG files into DIR instead of next to the SVG files.
const fs = require('fs'), path = require('path');
const {chromium} = require('playwright');
const args = process.argv.slice(2), oi = args.indexOf('--out');
const outDir = oi >= 0 ? args.splice(oi, 2)[1] : null;
(async () => {
  const browser = await chromium.launch();
  for (const f of args) {
    const svg = fs.readFileSync(f, 'utf8');
    const w = Math.ceil(+(svg.match(/<svg[^>]*\swidth="([\d.]+)/) || [])[1]), h = Math.ceil(+(svg.match(/<svg[^>]*\sheight="([\d.]+)/) || [])[1]);
    if (!w || !h) throw new Error(f + ': no width/height on <svg>');
    const page = await browser.newPage({viewport: {width: w, height: h}, deviceScaleFactor: 300 / 96});
    await page.setContent(`<!doctype html><html><body style="margin:0;background:#fff">${svg}</body></html>`);
    const out = path.join(outDir || path.dirname(f), path.basename(f, '.svg') + '.png');
    await page.screenshot({path: out, clip: {x: 0, y: 0, width: w, height: h}});
    await page.close();
    console.log(out, `${Math.round(w * 300 / 96)} x ${Math.round(h * 300 / 96)} px`);
  }
  await browser.close();
})();
