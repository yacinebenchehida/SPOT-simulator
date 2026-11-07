// Loss of neutral genetic diversity with migration (Supporting Information S3.3): SPOT's neutral simulations against
// Wright's recursion for a continent-island model, written out here independently of sim.js:
//   each generation, F <- 1/(2 Ne) + (1 - 1/(2 Ne)) F (haploids: 1/N), then F <- (1 - m)^2 F; H_g/H_0 = 1 - F_g.
// Simulated value: mean of p(1 - p) over 100,000 neutral trajectories relative to p0(1 - p0), every generation, with
// migrants carrying the starting frequency p0 (the default source frequency). Also reports the value given by SPOT
// (sim.js hetRetained, shown on the page). Default sampling design (1935, 1955, 1975, 2025: 90 generations).
// Writes paper/validation/migration_diversity.tsv. Run from the repository root: node paper/validation/migration_diversity_check.js
const S = require('../../sim.js');
const fs = require('fs');
const BASE = {system: 'diploid', neByBlock: false, neBlocks: null, genTime: 1, Nf: 1000, Nm: 1000, noDrift: false,
  p0: 0.3, s: 0, h: 0.5, onset: 1935, years: [1935, 1955, 1975, 2025], n: [0, 0, 0, 0], perfect: true,
  depthEarly: 10, depthLast: 15, err: 0.005, minDP: 3, alpha: 0.05, m: 1, sMale: 1, decline: 1, mig: 0, pMig: -1,
  kSib: 1, maleFrac: 0, reps: 1000, seed: 7};
const SC = [];
for (const m of [0.001, 0.01, 0.05]) SC.push([`diploid, Nf = Nm = 50, m = ${m}`, {system: 'diploid', Nf: 50, Nm: 50, mig: m}]);
for (const m of [0.01, 0.05]) SC.push([`haploid, N = 100, m = ${m}`, {system: 'haploid', Nf: 100, Nm: 100, mig: m}]);
for (const m of [0.01, 0.05]) SC.push([`haplodiploid, Nf = Nm = 100, m = ${m}`, {system: 'haplodiploid', Nf: 100, Nm: 100, mig: m}]);
SC.push(['diploid, Ne per interval 1,000 / 50 / 1,000, m = 0.01', {system: 'diploid', Nf: 500, Nm: 500, neByBlock: true, neBlocks: [1000, 50, 1000], mig: 0.01}]);
const ne = (sys, Nf, Nm) => sys === 'haploid' ? Nf : sys === 'haplodiploid' ? 9 * Nf * Nm / (4 * Nm + 2 * Nf) : 4 * Nf * Nm / (Nf + Nm);
const lines = ['scenario\tm\tgeneration\tH_ratio_simulated\tse\tH_ratio_theory\tH_ratio_without_migration'];
console.log('scenario\tH90/H0 simulated (SE)\ttheory\tSPOT value\twithout migration\tstandardised difference at 90\tlargest |difference| over generations');
let worst = 0;
for (const [lab, o] of SC) {
  const P = S.prep(Object.assign({}, BASE, o)), R = S.makeRng(31), reps = 100000, gens = S.gensOf(P), G = gens[gens.length - 1];
  const h = new Float64Array(G + 1), h2 = new Float64Array(G + 1);
  for (let r = 0; r < reps; r++) {
    const tr = S.simulateRep(P, 0, R, true);
    for (let g = 0; g <= G; g++) {
      const pf = tr.traj[g], pm = tr.trajM[g];
      const p = o.system === 'haploid' ? pf : o.system === 'diploid' ? (pf + pm) / 2 : (2 * pf + pm) / 3;
      h[g] += p * (1 - p); h2[g] += (p * (1 - p)) ** 2;
    }
  }
  const H0 = P.p0 * (1 - P.p0);
  let F = 0, F0 = 0, maxd = 0, last;
  for (let g = 0; g <= G; g++) {
    const sim = h[g] / reps / H0, se = Math.sqrt((h2[g] / reps - (h[g] / reps) ** 2) / reps) / H0;
    lines.push([lab, o.mig, g, sim.toFixed(5), se.toFixed(5), (1 - F).toFixed(5), (1 - F0).toFixed(5)].join('\t'));
    maxd = Math.max(maxd, Math.abs(sim - (1 - F))); last = [sim, se, 1 - F, 1 - F0];
    // next generation: drift with this generation's breeding numbers, then migration
    let f = 1;
    if (o.neByBlock) { let j = 0; while (j < gens.length - 2 && g >= gens[j + 1]) j++; f = o.neBlocks[j] / ne(o.system, o.Nf, o.Nm); }
    const e = ne(o.system, Math.max(2, Math.round(o.Nf * f)), Math.max(2, Math.round(o.Nm * f))), d = o.system === 'haploid' ? 1 / e : 1 / (2 * e);
    F = (1 - o.mig) ** 2 * (d + (1 - d) * F); F0 = d + (1 - d) * F0;
  }
  worst = Math.max(worst, maxd);
  console.log([lab, `${last[0].toFixed(4)} (${last[1].toFixed(4)})`, last[2].toFixed(4), S.hetRetained(Object.assign({}, BASE, o)).toFixed(4),
    last[3].toFixed(4), ((last[0] - last[2]) / last[1]).toFixed(2), maxd.toFixed(4)].join('\t'));
}
console.log('largest |simulated - theory| over all scenarios and generations:', worst.toFixed(4));
fs.writeFileSync('paper/validation/migration_diversity.tsv', lines.join('\n') + '\n');
