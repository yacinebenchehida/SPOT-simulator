// False-positive rate of the drift-aware threshold WITH rescaling.
// For each design, SPOT sets its threshold from rescaled neutral simulations (20 independent runs, default 1,000 + 2,000
// simulations each); the realised false-positive rate of each threshold is measured on 50,000 independent neutral
// simulations run WITHOUT rescaling (every generation simulated). The sampling-only threshold z(1 - alpha) is
// evaluated on the same set. Four sampling dates equally spaced over G generations, 20/20/20/40 individuals,
// default sequencing settings, p0 = 0.1, breeding females = breeding males = 10,000.
// Run from the repository root: node paper/validation/rescale_fpr.js > paper/validation/rescale_fpr.tsv
const S = require('../../sim.js');
const base = {genTime: 1, Nf: 1e4, Nm: 1e4, noDrift: false, p0: 0.1, s: 0, h: 0.5, onset: 0, n: [20, 20, 20, 40], perfect: false,
  depthEarly: 10, depthLast: 15, err: 0.005, minDP: 3, alpha: 0.05, m: 1, sMale: 1, decline: 1, mig: 0, pMig: -1, kSib: 1,
  maleFrac: 0, reps: 1000, seed: 1};
const N_THRESHOLDS = 20, N_TEST = 50000;
function trendZ(x, n, t, idx) { let N = 0, X = 0, T = 0; for (const i of idx) { N += n[i]; X += x[i]; T += n[i] * t[i]; } if (N === 0) return 0;
  const pb = X / N, tb = T / N; if (pb <= 0 || pb >= 1) return 0; let a = 0, b = 0;
  for (const i of idx) { a += x[i] * (t[i] - tb); b += n[i] * (t[i] - tb) ** 2; } return b <= 0 ? 0 : a / Math.sqrt(pb * (1 - pb) * b); }
console.log(['G', 'system', 'Ne', 'lambda', 'fpr_mean', 'fpr_sd', 'fpr_min', 'fpr_max', 'fpr_sampling_only', 'threshold_mean', 'n_test'].join('\t'));
const t0 = Date.now();
for (const G of [5000, 20000, 60000]) for (const system of ['haploid', 'diploid', 'haplodiploid']) {
  const P = Object.assign({}, base, {system, years: [0, Math.round(G / 3), Math.round(2 * G / 3), G]});
  const lambda = S.rescale(P).lambda;
  const crit = []; for (let i = 0; i < N_THRESHOLDS; i++) crit.push(S.runPower(Object.assign({}, P, {seed: 1000 + i})).trend);   // {crit, gamma}
  const Q = S.prep(Object.assign({}, P, {noRescale: true})), R = S.makeRng(99991), t = S.gensOf(Q), ALL = t.map((_, i) => i);
  const z = new Float64Array(N_TEST); for (let i = 0; i < N_TEST; i++) { const o = S.simulateRep(Q, 0, R, false); z[i] = trendZ(o.x, o.n, t, ALL); }
  const f = crit.map(c => { let k = 0; for (const v of z) k += S.exceed(v, c.crit, c.gamma); return k / N_TEST; });   // ties count with weight gamma
  const m = f.reduce((a, b) => a + b, 0) / f.length, sd = Math.sqrt(f.reduce((a, b) => a + (b - m) ** 2, 0) / (f.length - 1));
  const zc = S.qnorm(0.95); let naive = 0; for (const v of z) if (v > zc) naive++;
  console.log([G, system, Math.round(S.neFor(P)), lambda, m.toFixed(4), sd.toFixed(4), Math.min(...f).toFixed(4), Math.max(...f).toFixed(4),
    (naive / N_TEST).toFixed(4), (crit.reduce((a, b) => a + b.crit, 0) / crit.length).toFixed(3), N_TEST].join('\t'));
  console.error(G, system, 'lambda', lambda, 'FPR', m.toFixed(4), ((Date.now() - t0) / 1000).toFixed(0) + 's');
}
