// SLiM comparison, v2 (round numbers; repeated power comparison). Run from the repository root after
// paper/validation/slim_runs.sh and paper/validation/slim_runs_v2.sh:   node paper/validation/slim_compare_v2.js
// Writes:
//   paper/validation/slim_power_runs.tsv : power of each program in each of 20 independent runs per scenario with selection.
//        Run k = 1 uses the original SLiM runs (slim_out/<name>.tsv); runs k = 2..20 use slim_out/power/<name>_r<k>.tsv.
//        In each run and each program, the threshold is the 95% quantile of Z in that program's 2,000 neutral
//        simulations (same rule as SPOT), and power is the share of its 2,000 selection simulations above it.
//   paper/validation/slim_power_v2.tsv   : per scenario, mean and SD of power over the 20 runs for each program, the
//        difference of means (SPOT - SLiM) and its standard error from the between-run SDs.
//   paper/validation/slim_fpr_v2.tsv     : per neutral scenario, a threshold set from 20,000 SPOT neutral simulations and the
//        sampling-only threshold (1.645) applied to the 10,000 SLiM neutral simulations (original + _b, _c, _d, _e)
//        and to 10,000 independent SPOT neutral simulations.
// (The comparison of allele-frequency distributions is unchanged: paper/validation/slim_compare.js, slim_results.tsv.)
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

// ---- (2) power, 20 independent runs per program ----
const runs = [['scenario', 'run', 'threshold_spot', 'threshold_slim', 'power_spot', 'power_slim']];
const summ = [['scenario', 'system', 'NF', 'NM', 'p0', 's', 'h', 'runs', 'mean_power_spot', 'sd_power_spot',
  'mean_power_slim', 'sd_power_slim', 'difference', 'se_difference', 'z']];
SC.forEach(([sel, neu, system, Nf, Nm, p0, s, h], j) => {
  const ps = [], pl = [];
  for (let k = 1; k <= RUNS; k++) {
    const zSelSlim = Z(readSlim(slimFile(sel, k))), zNeuSlim = Z(readSlim(slimFile(neu, k)));
    const zSelSpot = Z(runSpot(system, Nf, Nm, p0, s, h, 10000 + 100 * j + k)), zNeuSpot = Z(runSpot(system, Nf, Nm, p0, 0, h, 20000 + 100 * j + k));
    if (zSelSlim.length !== R || zNeuSlim.length !== R) throw new Error(`${sel} run ${k}: wrong number of SLiM replicates`);
    const cS = quant95(zNeuSpot), cL = quant95(zNeuSlim), a = above(zSelSpot, cS), b = above(zSelSlim, cL);
    ps.push(a); pl.push(b);
    runs.push([sel, k, cS.crit.toFixed(3), cL.crit.toFixed(3), a.toFixed(4), b.toFixed(4)]);
  }
  const d = mean(ps) - mean(pl), se = Math.sqrt(sd(ps) ** 2 / RUNS + sd(pl) ** 2 / RUNS);
  summ.push([sel, system, Nf, Nm, p0, s, h, RUNS, mean(ps).toFixed(4), sd(ps).toFixed(4), mean(pl).toFixed(4), sd(pl).toFixed(4),
    d.toFixed(4), se.toFixed(4), (d / se).toFixed(2)]);
  console.error(sel, 'SPOT', mean(ps).toFixed(3), '+-', sd(ps).toFixed(3), 'SLiM', mean(pl).toFixed(3), '+-', sd(pl).toFixed(3));
});
fs.writeFileSync('paper/validation/slim_power_runs.tsv', runs.map(r => r.join('\t')).join('\n') + '\n');
fs.writeFileSync('paper/validation/slim_power_v2.tsv', summ.map(r => r.join('\t')).join('\n') + '\n');

// ---- (3) false-positive rates on 10,000 SLiM and 10,000 independent SPOT neutral simulations ----
const frow = [['neutral_scenario', 'system', 'NF', 'NM', 'p0', 'threshold_spot_20000', 'n_slim', 'fpr_spot_threshold_on_slim',
  'se', 'fpr_spot_threshold_on_spot', 'fpr_sampling_only_on_slim', 'fpr_sampling_only_on_spot', 'sd_Z_slim', 'sd_Z_spot']];
const done = new Set(); let q = 0;
for (const [, neu, system, Nf, Nm, p0, , h] of SC) {
  if (done.has(neu)) continue; done.add(neu); q++;
  const zSlim = ['', '_b', '_c', '_d', '_e'].flatMap(suf => Z(readSlim(neu + suf)));
  if (zSlim.length !== 10000) throw new Error(`${neu}: ${zSlim.length} SLiM neutral replicates instead of 10,000`);
  const big = []; for (let i = 0; i < 10; i++) big.push(...Z(runSpot(system, Nf, Nm, p0, 0, h, 30000 + 100 * q + i)));
  const test = []; for (let i = 0; i < 5; i++) test.push(...Z(runSpot(system, Nf, Nm, p0, 0, h, 40000 + 100 * q + i)));
  const c = quant95(big), f = above(zSlim, c);
  frow.push([neu, system, Nf, Nm, p0, c.crit.toFixed(3), zSlim.length, f.toFixed(4), Math.sqrt(0.05 * 0.95 / zSlim.length).toFixed(4),
    above(test, c).toFixed(4), above(zSlim, ZS).toFixed(4), above(test, ZS).toFixed(4), sd(zSlim).toFixed(3), sd(test).toFixed(3)]);
  console.error(neu, 'FPR SPOT threshold on SLiM', f.toFixed(4), 'on SPOT', above(test, c).toFixed(4), 'sampling-only SLiM', above(zSlim, ZS).toFixed(4), 'SPOT', above(test, ZS).toFixed(4));
}
fs.writeFileSync('paper/validation/slim_fpr_v2.tsv', frow.map(r => r.join('\t')).join('\n') + '\n');
