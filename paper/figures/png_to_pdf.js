//////////////////////////////////////////////////////////////////////////////////////////////////////
// Make a PDF from the PNG of a figure, without loss of quality                                     //
// The page of the PDF has the size of the image at 300 dpi.                                        //
//                                                                                                  //
// Needs Playwright (Chromium).                                                                     //
// Run from the repository root (the last two numbers are the width and height of                   //
// the image in pixels):                                                                            //
//   node paper/figures/png_to_pdf.js paper/figures/Figure1.png paper/figures/Figure1.pdf 3270 2771 //
//////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////
// Load the dependencies //
///////////////////////////
const {chromium} = require('playwright');    // Browser used to make the PDF
const path = require('path');
const fs = require('fs');                    // To read the PNG file


////////////////////////////////////////////////////////////////////////
// Show the image in the browser at its real size and print it to PDF //
////////////////////////////////////////////////////////////////////////
(async () => {
  const [src, out, W, H] = process.argv.slice(2);   // PNG file, PDF file, width and height in pixels
  const win = W / 300,                       // Width in inches at 300 dpi
    hin = H / 300;                           // Height in inches at 300 dpi
  const b = await chromium.launch();         // Start the browser
  const p = await b.newPage();
  const data = fs.readFileSync(src).toString('base64');   // Content of the PNG, written as text to put it in the page
  await p.setContent(`<html><body style="margin:0"><img src="data:image/png;base64,${data}" style="width:${win}in;height:${hin}in;display:block"></body></html>`);   // Page with only the image
  await p.pdf({path: out, width: win + 'in', height: hin + 'in', printBackground: true, pageRanges: '1'});   // Save one page of the size of the image
  await b.close();                           // Stop the browser
})();
