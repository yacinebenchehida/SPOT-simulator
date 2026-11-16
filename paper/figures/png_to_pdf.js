// Figure as a PDF from its PNG, without loss (page = image at 300 dpi). Needs Playwright (Chromium).
// Run from the repository root: node paper/figures/png_to_pdf.js paper/figures/Figure1.png paper/figures/Figure1.pdf 3270 2771
const {chromium} = require('playwright'); const path = require('path'); const fs = require('fs');
(async () => { const [src, out, W, H] = process.argv.slice(2); const win = W / 300, hin = H / 300;
  const b = await chromium.launch(); const p = await b.newPage();
  const data = fs.readFileSync(src).toString('base64');
  await p.setContent(`<html><body style="margin:0"><img src="data:image/png;base64,${data}" style="width:${win}in;height:${hin}in;display:block"></body></html>`);
  await p.pdf({path: out, width: win + 'in', height: hin + 'in', printBackground: true, pageRanges: '1'}); await b.close(); })();
