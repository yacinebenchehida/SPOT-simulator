// Depth case study (museum insect reference case): the same design sequenced at very low vs very high depth,
// with higher per-read error (e.g. post-mortem damage) at low depth, and unequal depth between dates.
// Everything else is the reference case. 5,000 selection and 10,000 neutral simulations per design.
// Run from the repository root: node paper/validation/depth_case.js > paper/worked_examples/depth_case.tsv
const {S, with: W} = require('./base.js');
const NE = [[1000, 1500], [66667, 100000]];            // breeding females = males -> Ne (haplodiploid)
const SS = [0.015, 0.02, 0.03];
const D = [];                                            // [label, depthEarly, depthLast, err, perfect, minDP]
D.push(['perfect', 10, 10, 0.005, true, 1]);
for (const d of [1, 2, 5, 10, 30]) D.push([`${d}x all dates`, d, d, 0.005, false, 1]);
for (const e of [0.02, 0.05]) for (const d of [1, 2]) D.push([`${d}x all dates, error ${e}`, d, d, e, false, 1]);
D.push(['2x earlier, 30x last', 2, 30, 0.005, false, 1]);
D.push(['30x earlier, 2x last', 30, 2, 0.005, false, 1]);
for (const d of [1, 2]) D.push([`${d}x all dates, minimum depth 3`, d, d, 0.005, false, 3]);
const H = ['Ne','s','design','depth_earlier','depth_last','error','min_depth','perfect','threshold','power','power_sampling_only',
  'fpr_sampling_only','missing','het_called_correctly','haploid_called_correctly','n_selection','n_neutral','seed'];
console.log(H.join('\t'));
let seed = 1000;
for (const [Nf, Ne] of NE) for (const s of SS) for (const [lab, de, dl, err, perfect, minDP] of D) {
  seed++;
  const r = S.runPower(W({Nf, Nm: Nf, s, depthEarly: de, depthLast: dl, err, perfect, minDP, reps: 5000, seed}));
  console.log([Ne, s, lab, de, dl, err, minDP, perfect, r.trend.crit.toFixed(3), r.trend.powerEmp.toFixed(4),
    r.trend.powerNaive.toFixed(4), r.trend.fprNaive.toFixed(4), r.missRate.toFixed(4),
    isNaN(r.hetAcc) ? '' : r.hetAcc.toFixed(4), isNaN(r.hapAcc) ? '' : r.hapAcc.toFixed(4), 5000, r.nNeutral, seed].join('\t'));
}
