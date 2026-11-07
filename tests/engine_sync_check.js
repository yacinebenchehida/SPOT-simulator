// Checks that the engine of the web page (index.html, <script id="spot-core">) is exactly sim.js, so that every analysis
// run with sim.js gives the same numbers as the page. Run from the repository root: node tests/engine_sync_check.js
const fs = require('fs');
const html = fs.readFileSync('index.html', 'utf8'), sim = fs.readFileSync('sim.js', 'utf8');
const tag = '<script id="spot-core">', a = html.indexOf(tag) + tag.length, core = html.slice(a, html.indexOf('</script>', a));
if (core.trim() === sim.trim()) { console.log('OK: the engine of index.html is identical to sim.js'); }
else { console.log('DIFFERENT: the engine of index.html is not sim.js'); process.exit(1); }
