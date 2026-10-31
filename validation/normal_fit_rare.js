// Normal-fit threshold beyond 1,000 loci (Supplement): with more than 1,000 candidate loci at alpha = 0.05, SPOT
// sets the threshold from a normal distribution fitted to 200,000 neutral simulations. This script measures its
// realised false-positive rate, for a rare (p0 = 0.01) and a common (p0 = 0.1) allele, in the three inheritance
// systems (Ne = 1,500; design of the reference case), on 2,000,000 independent neutral simulations, and compares
// power with that given by the empirical threshold of those 2,000,000 simulations.
// Run from the repository root: node validation/normal_fit_rare.js > validation/normal_fit_rare.tsv
const {S, with: W} = require('./base.js');
const SYS = [['haploid', 1500, 1500], ['diploid', 750, 750], ['haplodiploid', 1000, 1000]];
const NT = 2000000, MS = [2000, 10000];
console.log(['system','p0','m','alpha_per_locus','test','threshold_normal_fit','threshold_empirical_2e6','fpr_normal_fit',
  'ratio_to_alpha','power_normal_fit','power_empirical','n_neutral_fit','n_test'].join('\t'));
let seed = 3000;
for (const [system, Nf, Nm] of SYS) for (const p0 of [0.01, 0.1]) {
  const base = W({system, Nf, Nm, p0, s: 0.05});
  const prep = S.prep(S.rescale(base).P), t = S.gensOf(prep), K = t.length, ALL = t.map((_, i) => i), TWO = [0, K - 1];
  const Z = (x, n, idx) => { let N=0,X=0,T=0; for (const i of idx){N+=n[i];X+=x[i];T+=n[i]*t[i]} if(!N) return 0;
    const pb=X/N, tb=T/N; if(pb<=0||pb>=1) return 0; let a=0,b=0; for (const i of idx){a+=x[i]*(t[i]-tb); b+=n[i]*(t[i]-tb)**2}
    return b>0 ? a/Math.sqrt(pb*(1-pb)*b) : 0; };
  const R = S.makeRng(++seed), test = {trend: new Float64Array(NT), two: new Float64Array(NT)};
  for (let i = 0; i < NT; i++) { const o = S.simulateRep(prep, 0, R, false); test.trend[i] = Z(o.x, o.n, ALL); test.two[i] = Z(o.x, o.n, TWO); }
  const sorted = {trend: Float64Array.from(test.trend).sort(), two: Float64Array.from(test.two).sort()};
  for (const m of MS) {
    const r = S.runPower(Object.assign({}, base, {m, seed: ++seed})), a = 0.05 / m;
    for (const k of ['trend', 'two']) {
      const c = r[k].crit, ce = sorted[k][Math.ceil((1 - a) * NT) - 1];
      let f = 0; for (const z of test[k]) if (z > c) f++;
      let pe = 0; for (const z of r.sel[k]) if (z > ce) pe++;
      console.log([system, p0, m, a.toExponential(1), k, c.toFixed(3), ce.toFixed(3), (f / NT).toExponential(2), ((f / NT) / a).toFixed(1),
        r[k].powerEmp.toFixed(3), (pe / r.sel[k].length).toFixed(3), r.nNeutral, NT].join('\t'));
    }
  }
  console.error(system, p0, 'done');
}
