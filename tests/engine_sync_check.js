///////////////////////////////////////////////////////////////////////////////////
// Test: the simulation code inside the application is exactly sim.js            //
//                                                                               //
// The application (index.html) contains a copy of the simulation code, in the   //
// block <script id="spot-core">. This test checks that the copy is identical to //
// sim.js, so that every analysis run with sim.js gives the same numbers as the  //
// application.                                                                  //
//                                                                               //
// Run from the repository root:                                                 //
//   node tests/engine_sync_check.js                                             //
///////////////////////////////////////////////////////////////////////////////////


////////////////////////
// Read the two files //
////////////////////////
const fs = require('fs');
const html = fs.readFileSync('index.html', 'utf8'),   // The application
  sim = fs.readFileSync('sim.js', 'utf8');            // The simulation code


////////////////////////////////////////////////////
// Cut the simulation code out of the application //
////////////////////////////////////////////////////
const tag = '<script id="spot-core">',                // Start of the block in index.html
  a = html.indexOf(tag) + tag.length,                 // Position just after the start of the block
  core = html.slice(a, html.indexOf('</script>', a)); // Text of the block, up to its end


////////////////////////
// Compare and report //
////////////////////////
if (core.trim() === sim.trim()) {                     // Identical (spaces at both ends ignored)
  console.log('OK: the engine of index.html is identical to sim.js');
}
else {
  console.log('DIFFERENT: the engine of index.html is not sim.js');
  process.exit(1);                                    // Signal the failure to the calling program
}
