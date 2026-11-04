// No-drift example (Results 3.1): the two thresholds still differ slightly without drift, with unrelated individuals
// and perfect genotypes, when the favoured allele is rare (few distinct values of Z).
// Diploid, p0 = 0.01, s = 0.05, sampling design of the reference case. One SPOT run with the default numbers of
// simulations (1,000 selection, 2,000 neutral; seed 1) gives the threshold and the two powers; the false-positive
// rate of each threshold is then measured on 100,000 independent neutral simulations.
// Run from the repository root: node validation/nodrift_example.js > validation/nodrift_example.tsv
const {S, with: W} = require('./base.js');
const P = W({system: 'diploid', noDrift: true, perfect: true, p0: 0.01, s: 0.05, reps: 1000, seed: 1});
const r = S.runPower(P);
const crit = r.trend.crit, gamma = r.trend.gamma, ties = r.trend.ties, zc = r.zNaive;   // gamma: weight of ties at the threshold (sim.js)
const prep = S.prep(S.rescale(P).P), R = S.makeRng(777), t = S.gensOf(prep), all = t.map((_, i) => i);
const trendZ = (x, n) => { let N=0,X=0,T=0; for (const i of all){N+=n[i];X+=x[i];T+=n[i]*t[i]} if(!N) return 0;
  const pb=X/N, tb=T/N; if(pb<=0||pb>=1) return 0; let a=0,b=0; for (const i of all){a+=x[i]*(t[i]-tb); b+=n[i]*(t[i]-tb)**2}
  return b>0 ? a/Math.sqrt(pb*(1-pb)*b) : 0; };
// power of the neutral-simulation threshold with three ways of counting the selection simulations tied with it
let sGt = 0, sEq = 0; for (const z of r.sel.trend) { if (S.isTie(z, crit)) sEq++; else if (z > crit) sGt++; }
const NT = 100000; let f1 = 0, f2 = 0;
for (let i = 0; i < NT; i++) { const o = S.simulateRep(prep, 0, R, false); const z = trendZ(o.x, o.n); f1 += S.exceed(z, crit, gamma); if (z > zc) f2++; }
console.log(['threshold_neutral_simulations','n_ties_at_threshold','gamma_ties','threshold_sampling_only','power_neutral_simulation_threshold','power_sampling_only',
  'fpr_neutral_simulation_threshold','fpr_sampling_only','n_independent_neutral','n_selection_above','n_selection_tied','power_strictly_above','power_ties_all_counted'].join('\t'));
console.log([crit.toFixed(3), ties, gamma.toFixed(4), zc.toFixed(3), r.trend.powerEmp.toFixed(4), r.trend.powerNaive.toFixed(4), (f1/NT).toFixed(4), (f2/NT).toFixed(4), NT, sGt, sEq, (sGt/P.reps).toFixed(4), ((sGt+sEq)/P.reps).toFixed(4)].join('\t'));
