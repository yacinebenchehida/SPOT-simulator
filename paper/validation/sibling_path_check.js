// Check of the family (sibling) sampling procedure: with families of one individual (k = 1), sampling through the family
// branch of sim.js (P.famPath) must give the same power as the default sampling of unrelated individuals.
// In the default procedure, each individual is drawn from the genotype frequencies of the offspring after selection; in
// the family procedure, a mother and a father are drawn from the female and male offspring after selection, and the
// individual is their offspring (Mendelian inheritance, no further selection). Both procedures give unrelated individuals.
// Diploid and haplodiploid, Ne = 300, 1,500 and 100,000, s = 0.02 and 0.05, design of the reference case (p0 = 0.1,
// 4 sampling dates, 20/20/20/40 females), perfect genotypes, 20,000 selection and 40,000 neutral simulations per run,
// 3 seeds. Each run sets its own threshold, as SPOT does.
// se_difference: binomial standard error of the difference between two powers from 20,000 selection simulations each
// (it ignores the variation of the threshold between runs, so it slightly underestimates the true SE).
// Run from the repository root: node paper/validation/sibling_path_check.js > paper/validation/sibling_path_check.tsv  (about 5 minutes)
const {S, with: W} = require('./base.js');
const R = 20000, SEEDS = [1, 2, 3];
const sizes = (system, Ne) => system === 'diploid' ? {Nf: Ne / 2, Nm: Ne / 2} : {Nf: Ne / 1.5, Nm: Ne / 1.5};
console.log(['system', 'Ne', 's', 'seed', 'statistic', 'power_default', 'power_family_path', 'difference', 'se_difference', 'z'].join('\t'));
const t0 = Date.now();
for (const system of ['diploid', 'haplodiploid']) for (const Ne of [300, 1500, 100000]) for (const s of [0.02, 0.05]) for (const seed of SEEDS) {
  const P = W(Object.assign({system, s, perfect: true, reps: R, seed}, sizes(system, Ne)));
  const a = S.runPower(P), b = S.runPower(Object.assign({}, P, {famPath: true}));
  for (const k of ['trend', 'two']) {
    const p1 = a[k].powerEmp, p2 = b[k].powerEmp, d = p2 - p1, se = Math.sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / R);
    console.log([system, Ne, s, seed, k, p1.toFixed(4), p2.toFixed(4), d.toFixed(4), se.toFixed(4), se > 0 ? (d / se).toFixed(2) : ''].join('\t'));
  }
  console.error(system, Ne, s, seed, ((Date.now() - t0) / 1000).toFixed(0) + 's');
}
