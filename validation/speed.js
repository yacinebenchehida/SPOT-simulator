// Speed of one complete calculation (simulations, thresholds and power of both tests), as in Supporting Information S5:
// the five example scenarios of the page and the reference case, with the default numbers of simulations
// (1,000 selection, 2,000 neutral) and one locus, and the reference case with 1,000 candidate loci (200,000 neutral).
// Each calculation is timed 5 times (seeds 1-5) after one warm-up run; one processor core (Node.js runs single-threaded);
// run with no other calculation running on the computer.
// Run from the repository root: node validation/speed.js > validation/speed.tsv
const os = require('os');
const S = require('../sim.js');
const {REF} = require('./base.js');
const D = {system: 'haplodiploid', neByBlock: false, neBlocks: null, genTime: 1, Nf: 1000, Nm: 1000, noDrift: false,
  p0: 0.01, s: 0.05, h: 0.5, onset: 1935, years: [1935, 1955, 1975, 2025], n: [20, 20, 20, 40], perfect: false,
  depthEarly: 10, depthLast: 15, err: 0.005, minDP: 3, alpha: 0.05, m: 1, sMale: 1, decline: 1, mig: 0, pMig: -1,
  kSib: 1, maleFrac: 0, reps: 1000, seed: 1};
// the example scenarios of the page (same values as in monte_carlo_precision.js)
const SC = [
  ['Bumblebees', D],
  ['Long-lived mammal', Object.assign({}, D, {system: 'diploid', genTime: 20, Nf: 2500, Nm: 2500, neByBlock: true, neBlocks: [5000, 200], p0: 0.1, s: 0.2, onset: 1880, years: [1880, 1950, 2020], n: [20, 20, 40], depthEarly: 4, depthLast: 15})],
  ['Human ancient DNA', Object.assign({}, D, {system: 'diploid', genTime: 28, Nf: 5000, Nm: 5000, p0: 0.1, s: 0.005, onset: -7000, years: [-7000, -5000, -3000, -1000, 1000, 2000], n: [50, 50, 50, 50, 50, 100], depthEarly: 2, depthLast: 20, err: 0.01, minDP: 1})],
  ['Bacterial pathogen', Object.assign({}, D, {system: 'haploid', genTime: 0.01, Nf: 100000, Nm: 100000, p0: 0.05, s: 0.0002, onset: 1850, years: [1850, 1900, 1950, 2000, 2020], n: [10, 10, 15, 20, 40], depthEarly: 10, depthLast: 30})],
  ['Large diploid population', Object.assign({}, D, {system: 'diploid', genTime: 1, Nf: 1e6, Nm: 1e6, p0: 0.05, s: 0.02, onset: 1900, years: [1900, 1950, 2000, 2020], n: [30, 30, 30, 60]})],
  ['Reference case', REF],
  ['Reference case, 1,000 candidate loci', Object.assign({}, REF, {m: 1000})],
];
const cpu = os.cpus()[0].model.trim(), node = process.version;
console.log(['scenario', 'loci', 'selection_simulations', 'neutral_simulations', 'median_s', 'min_s', 'max_s', 'cpu', 'node'].join('\t'));
for (const [lab, P] of SC) {
  S.runPower(Object.assign({}, P, {seed: 99}));            // warm-up
  const t = [];
  for (let k = 1; k <= 5; k++) { const a = process.hrtime.bigint(); S.runPower(Object.assign({}, P, {seed: k})); t.push(Number(process.hrtime.bigint() - a) / 1e9); }
  t.sort((x, y) => x - y);
  console.log([lab, P.m, P.reps, S.neutralReps(P), t[2].toFixed(3), t[0].toFixed(3), t[4].toFixed(3), cpu, node].join('\t'));
}
