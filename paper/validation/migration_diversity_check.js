////////////////////////////////////////////////////////////////////////////////////
// Loss of neutral genetic diversity with migration (Supporting Information S3)   //
//                                                                                //
// The neutral simulations of SPOT are compared with theory (Wright's recursion   //
// for a continent-island model), written here independently of sim.js:           //
//   each generation, F <- 1/(2 Ne) + (1 - 1/(2 Ne)) F   (haploids: 1/N)          //
//   then            F <- (1 - m)^2 F                                             //
//   and the diversity left is H_g / H_0 = 1 - F_g.                               //
//                                                                                //
// Simulated value: the mean of p(1 - p) over 100,000 neutral trajectories,       //
// relative to p0(1 - p0), at every generation. Migrants carry the starting       //
// frequency p0 (the default). The script also reports the value given by SPOT    //
// (hetRetained in sim.js, shown in the application).                             //
// Sampling dates of the reference case (1935, 1955, 1975, 2025: 90 generations). //
//                                                                                //
// Run from the repository root:                                                  //
//   node paper/validation/migration_diversity_check.js                           //
// It writes paper/validation/migration_diversity.tsv.                            //
////////////////////////////////////////////////////////////////////////////////////


///////////////////////////
// Load the dependencies //
///////////////////////////
const S = require('../../sim.js');           // Simulation code of SPOT
const fs = require('fs');                    // To write the results file


//////////////////////////////////////////////////////////////////
// Settings shared by all scenarios (no selection, no sampling) //
//////////////////////////////////////////////////////////////////
const BASE = {system: 'diploid',             // Inheritance system (changed below)
  neByBlock: false, neBlocks: null,          // Ne is not set separately for each interval (changed in the last scenario)
  genTime: 1,                                // Years per generation
  Nf: 1000, Nm: 1000,                        // Breeding females and males (changed below)
  noDrift: false,                            // Drift is simulated
  p0: 0.3,                                   // Starting frequency of the allele
  s: 0,                                      // No selection
  h: 0.5,                                    // Dominance
  onset: 1935,                               // Year selection starts (not used, s = 0)
  years: [1935, 1955, 1975, 2025],           // Sampling years (90 generations)
  n: [0, 0, 0, 0],                           // No individual sampled: only the population matters here
  perfect: true,                             // Perfect genotypes
  depthEarly: 10, depthLast: 15, err: 0.005, minDP: 3,   // Sequencing settings (not used)
  alpha: 0.05, m: 1,                         // Significance level and number of loci (not used)
  sMale: 1, decline: 1,                      // No sex difference in selection, constant Ne
  mig: 0, pMig: -1,                          // Migration rate (changed below); migrants carry the starting frequency
  kSib: 1, maleFrac: 0,                      // Unrelated individuals, females only
  reps: 1000,                                // Not used here
  seed: 7};                                  // Random seed


////////////////////////////////////////////////////////////
// The scenarios: [label, settings that differ from BASE] //
////////////////////////////////////////////////////////////
const SC = [];
for (const m of [0.001, 0.01, 0.05])         // Diploid, three migration rates
  SC.push([`diploid, Nf = Nm = 50, m = ${m}`, {system: 'diploid', Nf: 50, Nm: 50, mig: m}]);
for (const m of [0.01, 0.05])                // Haploid, two migration rates
  SC.push([`haploid, N = 100, m = ${m}`, {system: 'haploid', Nf: 100, Nm: 100, mig: m}]);
for (const m of [0.01, 0.05])                // Haplodiploid, two migration rates
  SC.push([`haplodiploid, Nf = Nm = 100, m = ${m}`, {system: 'haplodiploid', Nf: 100, Nm: 100, mig: m}]);
SC.push(['diploid, Ne per interval 1,000 / 50 / 1,000, m = 0.01', {system: 'diploid', Nf: 500, Nm: 500, neByBlock: true, neBlocks: [1000, 50, 1000], mig: 0.01}]);   // A bottleneck with migration


///////////////////////////////////////////////////////////////////
// Function: effective population size from the breeding numbers //
// (written here independently of sim.js)                        //
///////////////////////////////////////////////////////////////////
const ne = (sys, Nf, Nm) => sys === 'haploid' ? Nf                    // Haploid: N
  : sys === 'haplodiploid' ? 9 * Nf * Nm / (4 * Nm + 2 * Nf)          // Haplodiploid
  : 4 * Nf * Nm / (Nf + Nm);                                          // Diploid


////////////////////////////////
// Headers of the two outputs //
////////////////////////////////
const lines = ['scenario\tm\tgeneration\tH_ratio_simulated\tse\tH_ratio_theory\tH_ratio_without_migration'];   // Results file: one line per scenario and generation
console.log('scenario\tH90/H0 simulated (SE)\ttheory\tSPOT value\twithout migration\tstandardised difference at 90\tlargest |difference| over generations');   // Summary on the screen


////////////////////////////////////////////////////////////////////
// For each scenario: simulate, then compare with theory at every //
// generation                                                     //
////////////////////////////////////////////////////////////////////
let worst = 0;                               // Largest difference between simulation and theory, over all scenarios
for (const [lab, o] of SC) {
  const P = S.prep(Object.assign({}, BASE, o)),   // Settings of this scenario, in the form used by the simulation code
    R = S.makeRng(31),                       // Random number generator
    reps = 100000,                           // Number of neutral trajectories
    gens = S.gensOf(P),                      // Generation number of each sampling date
    G = gens[gens.length - 1];               // Last generation (90)

  // Simulate the trajectories and add up p(1 - p) at every generation
  const h = new Float64Array(G + 1),         // Sum of p(1 - p) at each generation
    h2 = new Float64Array(G + 1);            // Sum of its square (for the standard error)
  for (let r = 0; r < reps; r++) {
    const tr = S.simulateRep(P, 0, R, true); // One neutral trajectory (allele frequency in females and in males)
    for (let g = 0; g <= G; g++) {
      const pf = tr.traj[g],                 // Frequency in females at generation g
        pm = tr.trajM[g];                    // Frequency in males at generation g
      const p = o.system === 'haploid' ? pf : o.system === 'diploid' ? (pf + pm) / 2 : (2 * pf + pm) / 3;   // Frequency in the population (haplodiploid: females carry two thirds of the gene copies)
      h[g] += p * (1 - p);
      h2[g] += (p * (1 - p)) ** 2;
    }
  }

  // Compare with Wright's recursion, generation by generation
  const H0 = P.p0 * (1 - P.p0);              // Diversity at the first date
  let F = 0,                                 // Theory with migration
    F0 = 0,                                  // Theory without migration
    maxd = 0,                                // Largest difference between simulation and theory in this scenario
    last;                                    // Values at the last generation
  for (let g = 0; g <= G; g++) {
    const sim = h[g] / reps / H0,            // Simulated diversity left at generation g
      se = Math.sqrt((h2[g] / reps - (h[g] / reps) ** 2) / reps) / H0;   // Its standard error
    lines.push([lab, o.mig, g, sim.toFixed(5), se.toFixed(5), (1 - F).toFixed(5), (1 - F0).toFixed(5)].join('\t'));
    maxd = Math.max(maxd, Math.abs(sim - (1 - F)));
    last = [sim, se, 1 - F, 1 - F0];
    // next generation: drift with this generation's breeding numbers, then migration
    let f = 1;                               // Ne of this generation relative to the Ne given by Nf and Nm
    if (o.neByBlock) {                       // Ne set separately for each interval between sampling dates
      let j = 0;
      while (j < gens.length - 2 && g >= gens[j + 1]) j++;   // Interval that contains generation g
      f = o.neBlocks[j] / ne(o.system, o.Nf, o.Nm);
    }
    const e = ne(o.system, Math.max(2, Math.round(o.Nf * f)), Math.max(2, Math.round(o.Nm * f))),   // Ne of this generation
      d = o.system === 'haploid' ? 1 / e : 1 / (2 * e);   // Diversity lost by drift in one generation
    F = (1 - o.mig) ** 2 * (d + (1 - d) * F);             // Drift, then migration
    F0 = d + (1 - d) * F0;                                // Drift only
  }
  worst = Math.max(worst, maxd);

  // Summary line of this scenario
  console.log([lab, `${last[0].toFixed(4)} (${last[1].toFixed(4)})`, last[2].toFixed(4), S.hetRetained(Object.assign({}, BASE, o)).toFixed(4),
    last[3].toFixed(4), ((last[0] - last[2]) / last[1]).toFixed(2), maxd.toFixed(4)].join('\t'));
}


//////////////////////////////////////
// Overall summary and results file //
//////////////////////////////////////
console.log('largest |simulated - theory| over all scenarios and generations:', worst.toFixed(4));
fs.writeFileSync('paper/validation/migration_diversity.tsv', lines.join('\n') + '\n');
