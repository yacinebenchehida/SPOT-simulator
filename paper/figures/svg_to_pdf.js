///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Make a vector PDF from the SVG of a figure                                                            //
// The page of the PDF has the size of the figure.                                                       //
//                                                                                                       //
// Needs Playwright (Chromium).                                                                          //
// Run from the repository root:                                                                         //
//   node paper/figures/svg_to_pdf.js paper/validation/false_positive_rate.svg paper/figures/Figure2.pdf //
///////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////
// Load the dependencies //
///////////////////////////
const {chromium} = require('playwright');    // Browser used to make the PDF
const fs = require('fs');                    // To read the SVG file


/////////////////////////////////////////////////////
// Show the SVG in the browser and print it to PDF //
/////////////////////////////////////////////////////
(async () => {
  const [src, out] = process.argv.slice(2);  // SVG file, PDF file
  const svg = fs.readFileSync(src, 'utf8');  // Content of the SVG file
  const [, W, H] = svg.match(/width="([\d.]+)" height="([\d.]+)"/);   // Width and height of the figure, read from the <svg> tag
  const b = await chromium.launch();         // Start the browser
  const p = await b.newPage();
  await p.setContent(`<html><body style="margin:0">${svg}</body></html>`);   // Page with only the figure
  await p.pdf({path: out, width: W + 'px', height: H + 'px', printBackground: true, pageRanges: '1'});   // Save one page of the size of the figure
  await b.close();                           // Stop the browser
})();
