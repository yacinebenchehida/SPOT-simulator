//////////////////////////////////////////////////////////////////////////////////////////
// Convert SVG figures to PNG at 300 dpi                                                //
//                                                                                      //
// An SVG pixel is 1/96 inch, so the image is enlarged by 300/96. The conversion        //
// uses the Chromium browser installed with Playwright.                                 //
//                                                                                      //
// Run from the repository root:                                                        //
//   node paper/validation/svg_to_png.js file1.svg [file2.svg ...]   ->   file1.png ... //
// Optional: --out DIR writes the PNG files into DIR instead of next to the SVG files.  //
//////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////
// Load the dependencies //
///////////////////////////
const fs = require('fs'),                    // To read the SVG files
  path = require('path');                    // To build the names of the PNG files
const {chromium} = require('playwright');    // Browser used to draw the SVG


///////////////////////////////////
// Arguments of the command line //
///////////////////////////////////
const args = process.argv.slice(2),          // The SVG files (and the optional --out DIR)
  oi = args.indexOf('--out');                // Position of --out, or -1 if absent
const outDir = oi >= 0 ? args.splice(oi, 2)[1] : null;   // Output folder; --out and DIR are removed from the list of files


///////////////////////////////////////////////////////////
// Draw each SVG in the browser and save a picture of it //
///////////////////////////////////////////////////////////
(async () => {
  const browser = await chromium.launch();   // Start the browser
  for (const f of args) {                    // Loop over the SVG files
    const svg = fs.readFileSync(f, 'utf8');  // Content of the SVG file
    const w = Math.ceil(+(svg.match(/<svg[^>]*\swidth="([\d.]+)/) || [])[1]),    // Width of the figure, read from the <svg> tag
      h = Math.ceil(+(svg.match(/<svg[^>]*\sheight="([\d.]+)/) || [])[1]);       // Height of the figure
    if (!w || !h) throw new Error(f + ': no width/height on <svg>');             // Stop if the size is missing
    const page = await browser.newPage({viewport: {width: w, height: h}, deviceScaleFactor: 300 / 96});   // Page of the size of the figure, enlarged to 300 dpi
    await page.setContent(`<!doctype html><html><body style="margin:0;background:#fff">${svg}</body></html>`);   // Show the SVG on a white background
    const out = path.join(outDir || path.dirname(f), path.basename(f, '.svg') + '.png');   // Name of the PNG file
    await page.screenshot({path: out, clip: {x: 0, y: 0, width: w, height: h}});            // Save the picture
    await page.close();
    console.log(out, `${Math.round(w * 300 / 96)} x ${Math.round(h * 300 / 96)} px`);       // Print the name and size of the PNG
  }
  await browser.close();                     // Stop the browser
})();
