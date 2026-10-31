// Reference case (museum insect), as in the manuscript (Table 3) — shared by the new analyses.
const S = require('../sim.js');
const REF = {system:'haplodiploid', neByBlock:false, neBlocks:null, genTime:1, Nf:1000, Nm:1000, noDrift:false,
  p0:0.1, s:0.02, h:0.5, onset:1935, years:[1935,1955,1975,2025], n:[20,20,20,40], perfect:false,
  depthEarly:10, depthLast:15, err:0.005, minDP:3, alpha:0.05, m:1, sMale:1, decline:1, mig:0, pMig:-1,
  kSib:1, maleFrac:0, reps:1000, seed:1};
module.exports = {S, REF, with: o => Object.assign({}, REF, o)};
