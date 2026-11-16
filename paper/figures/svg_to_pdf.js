// Figure as a vector PDF from its SVG (page = figure size). Needs Playwright (Chromium).
// Run from the repository root: node paper/figures/svg_to_pdf.js paper/validation/false_positive_rate.svg paper/figures/Figure2.pdf
const {chromium} = require('playwright'); const fs = require('fs');
(async () => { const [src, out] = process.argv.slice(2); const svg = fs.readFileSync(src, 'utf8');
  const [, W, H] = svg.match(/width="([\d.]+)" height="([\d.]+)"/);
  const b = await chromium.launch(); const p = await b.newPage();
  await p.setContent(`<html><body style="margin:0">${svg}</body></html>`);
  await p.pdf({path: out, width: W + 'px', height: H + 'px', printBackground: true, pageRanges: '1'}); await b.close(); })();
