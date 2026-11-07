// SLiM comparison, part 4: power that SPOT's threshold gives on SLiM's populations with selection.
// For each scenario with selection, a threshold is set from 20,000 SPOT neutral simulations (95% quantile, as in SPOT)
// and applied to the selection simulations of both programs in each of the 20 independent runs (2,000 per run):
// if SPOT simulates selection and drift as SLiM does, SPOT's threshold gives the same power on both.
// Run from the repository root after slim_runs.sh and slim_runs_v2.sh: node paper/validation/slim_power_spot_threshold.js
const S = require('../../sim.js');
const fs = require('fs');
const GENS = [0, 20, 40, 90], NS = [20, 20, 20, 40], R = 2000, RUNS = 20, ZS = 1.6448536269514722;
const SC = [
  ['diploid_N1000_s0.02', 'diploid_N1000_neutral', 'diploid', 500, 500, 0.1, 0.02, 0.5],
  ['diploid_N1000_s0.05_rec', 'diploid_N1000_neutral', 'diploid', 500, 500, 0.1, 0.05, 0],
  ['diploid_N200_s0.05', 'diploid_N200_neutral', 'diploid', 100, 100, 0.1, 0.05, 0.5],
  ['haploid_N1000_s0.02', 'haploid_N1000_neutral', 'haploid', 1000, 1000, 0.05, 0.02, 0.5],
  ['haplodiploid_N1000_s0.02', 'haplodiploid_N1000_neutral', 'haplodiploid', 500, 500, 0.1, 0.02, 0.5],
  ['haplodiploid_N2000_s0.05', 'haplodiploid_N2000_neutral', 'haplodiploid', 1000, 1000, 0.01, 0.05, 0.5],
];
// ---- reading SLiM output: columns replicate, generation, pfs, x, nall ----
function readSlim(name) {
  const reps = new Map();
  for (const line of fs.readFileSync(`paper/validation/slim_out/${name}.tsv`, 'utf8').trim().split('\n')) {
    const [id, g, pfs, x, n] = line.split('\t').map(Number), k = GENS.indexOf(g);
    if (!reps.has(id)) reps.set(id, {pfs: [], x: [], n: []});
    const r = reps.get(id); r.pfs[k] = pfs; r.x[k] = x; r.n[k] = n;
  }
  return [...reps.values()];
}

// ---- SPOT: same scenario with its own simulation code ----
function runSpot(system, Nf, Nm, p0, s, h, seed) {
  const P = S.prep({system, neByBlock: false, neBlocks: null, genTime: 1, Nf, Nm, noDrift: false, p0, s, h, onset: 0,
    years: GENS, n: NS, perfect: true, depthEarly: 10, depthLast: 15, err: 0.005, minDP: 3, alpha: 0.05, m: 1, sMale: 1,
    decline: 1, mig: 0, pMig: -1, kSib: 1, maleFrac: 0, reps: R, seed});
  const rng = S.makeRng(seed), out = [];
  for (let i = 0; i < R; i++) {
    const t = S.simulateRep(P, s, rng, true);
    out.push({pfs: GENS.map(g => t.traj[g]), x: t.x, n: t.n});
  }
  return out;
}

// ---- statistics (identical for both simulators) ----
function trendZ(x, n) {    // Cochran-Armitage trend statistic, as trendZ() in sim.js
  let N = 0, X = 0, T = 0;
  for (let i = 0; i < x.length; i++) { N += n[i]; X += x[i]; T += n[i] * GENS[i]; }
  const pb = X / N, tb = T / N; if (pb <= 0 || pb >= 1) return 0;
  let num = 0, den = 0;
  for (let i = 0; i < x.length; i++) { num += x[i] * (GENS[i] - tb); den += n[i] * (GENS[i] - tb) ** 2; }
  return num / Math.sqrt(pb * (1 - pb) * den);
}
// 95% quantile with SPOT's rule for ties at the threshold (empiricalThreshold in sim.js): {crit, gamma, ties}
const quant95 = z => S.empiricalThreshold(z, 0.05);
const mean = a => a.reduce((p, q) => p + q, 0) / a.length;
const sd = a => { const m = mean(a); return Math.sqrt(a.reduce((p, q) => p + (q - m) ** 2, 0) / (a.length - 1)); };
function ks(a, b) {        // two-sample Kolmogorov-Smirnov statistic D and asymptotic p-value
  const A = [...a].sort((p, q) => p - q), B = [...b].sort((p, q) => p - q);
  let i = 0, j = 0, D = 0;
  while (i < A.length && j < B.length) {
    const v = Math.min(A[i], B[j]);
    while (i < A.length && A[i] <= v) i++;
    while (j < B.length && B[j] <= v) j++;
    D = Math.max(D, Math.abs(i / A.length - j / B.length));
  }
  const ne = A.length * B.length / (A.length + B.length), lam = (Math.sqrt(ne) + 0.12 + 0.11 / Math.sqrt(ne)) * D;
  let p = 0; for (let k = 1; k <= 100; k++) p += 2 * (-1) ** (k - 1) * Math.exp(-2 * k * k * lam * lam);
  return [D, Math.min(1, Math.max(0, p))];
}


const Z = rs => rs.map(r => trendZ(r.x, r.n));
// share of Z above a threshold: SPOT threshold {crit, gamma} (ties count with weight gamma), or a fixed number (sampling-only)
const above = (z, c) => typeof c === 'number' ? mean(z.map(v => v > c ? 1 : 0)) : mean(z.map(v => S.exceed(v, c.crit, c.gamma)));
const slimFile = (name, k) => k === 1 ? name : `power/${name}_r${k}`;


const rows = [['scenario', 'threshold_spot_20000', 'runs', 'mean_power_on_spot', 'sd_power_on_spot', 'mean_power_on_slim', 'sd_power_on_slim', 'difference', 'se_difference', 'z']];
SC.forEach(([sel, neu, system, Nf, Nm, p0, s, h], j) => {
  const big = []; for (let i = 0; i < 10; i++) big.push(...Z(runSpot(system, Nf, Nm, p0, 0, h, 50000 + 100 * j + i)));
  const c = quant95(big), ps = [], pl = [];
  for (let k = 1; k <= RUNS; k++) {
    ps.push(above(Z(runSpot(system, Nf, Nm, p0, s, h, 10000 + 100 * j + k)), c));
    const zl = Z(readSlim(slimFile(sel, k))); if (zl.length !== R) throw new Error(sel + ' run ' + k);
    pl.push(above(zl, c));
  }
  const d = mean(ps) - mean(pl), se = Math.sqrt(sd(ps) ** 2 / RUNS + sd(pl) ** 2 / RUNS);
  rows.push([sel, c.crit.toFixed(3), RUNS, mean(ps).toFixed(4), sd(ps).toFixed(4), mean(pl).toFixed(4), sd(pl).toFixed(4), d.toFixed(4), se.toFixed(4), (d / se).toFixed(2)]);
  console.error(sel, 'SPOT threshold', c.crit.toFixed(3), 'power on SPOT', mean(ps).toFixed(3), 'on SLiM', mean(pl).toFixed(3), 'z', (d / se).toFixed(2));
});
require('fs').writeFileSync('paper/validation/slim_power_spot_threshold.tsv', rows.map(r => r.join('\t')).join('\n') + '\n');
