# Ties at the neutral-simulation threshold: what changed

## The change

When the score Z takes few distinct values (rare allele, little drift, small samples, first-versus-last test), several neutral simulations can have exactly the value of the threshold c. SPOT counted only values strictly above c, so the realised false-positive rate was at most α′, and sometimes well below it (1.7% instead of 5% in the no-drift example).

SPOT now shares the tied value deterministically (`empiricalThreshold()` and `exceed()` in `sim.js`):

- A = number of neutral Z > c, B = number of neutral Z equal to c (relative tolerance 1e-9 · max(1, |c|));
- γ = (α′·n − A)/B, limited to [0, 1); γ = 0 when B = 0;
- exactly a proportion α′ of the neutral simulations lies above the threshold: (A + γB)/n = α′;
- power = (#selection Z > c + γ · #selection Z = c) / number of selection simulations.

The threshold value itself, the normal-fit branch (more than 1,000 loci, or 100 for the power curve and planner) and the sampling-only threshold are unchanged. Every evaluation script counts the false positives of a SPOT threshold, and the power at a threshold, with the same rule.

## Files

| File | Contents |
|---|---|
| `ties_summary.tsv` | The numbers quoted in the manuscript, recomputed from the old (commit 37bb426) and the new result tables: old value, new value, script, where they appear. |
| `ties_changes.tsv` | Every cell that changed in every result table: file, script, row, column, old value, new value, difference (run times excluded). |

Old tables: `git show 37bb426:<file>`.

## Numbers quoted in the manuscript that change

| Quantity | Old | New | Script |
|---|---|---|---|
| False-positive rate, mean of 20 thresholds, trend test (36 cases) | 4.55–5.27% | 4.83–5.27% | `false_positive_rate.js` A, B |
| Same, first-versus-last test | 3.69–5.11% | 4.84–5.25% | `false_positive_rate.js` A, B |
| Single thresholds, first-versus-last test | 2.8–7.3% | 3.5–7.3% | `false_positive_rate.js` A, B |
| p₀ = 0.01, trend / first-versus-last (mean) | 4.55% / 4.28% | 4.92% / 5.00% | `false_positive_rate.js` B |
| Many loci, realised/intended, trend | 0.91–1.22 | 0.96–1.22 | `false_positive_rate.js` C |
| Many loci, realised/intended, first-versus-last | 0.67–1.16 | 0.98–1.18 | `false_positive_rate.js` C |
| SLiM: power with SPOT's threshold, SPOT − SLiM | −0.3 to +0.3 points | −0.4 to +0.3 points | `slim_power_spot_threshold.js` |
| SD of power, first-versus-last, 1,000 selection simulations | 1.5–3.1 points | 1.4–3.1 points | `monte_carlo_precision.js` |
| Rescaling (trend): higher / lower among 107 | 50 / 52 | 50 / 53 | `rescale_validation.js` |
| Rescaling (trend): SD of differences (without the outlier) | 2.3 (1.8) points | 2.4 (1.9) points | `rescale_validation.js` |
| Rescaling (trend): differences > 2 SE | 8 | 7 | `rescale_validation.js` |
| Allocation, Ne = 1,500, +60 historical (trend) | +2.0 to +4.6 points | +1.7 to +4.6 points | `worked_examples.js` B |
| Allocation, Ne = 1,500, +60 modern (trend) | +2.2 to +5.6 points | +2.1 to +5.6 points | `worked_examples.js` B |
| Allocation, Ne = 1,500, first and last only (trend) | −0.3 to +6.4 points | −0.5 to +6.5 points | `worked_examples.js` B |
| Allocation, Ne = 100,000, first and last only (trend) | +12.4 to +15.5 points | +12.4 to +15.7 points | `worked_examples.js` B |
| No-drift example: false positives / power | 1.65% / 76.80% | 4.85% / 77.02% | `nodrift_example.js` |

Other changes (Table S5 power, regression of SLiM on SPOT, first-versus-last columns of the case studies, normal-fit supplement, Table S7) are listed in `ties_summary.tsv` and `ties_changes.tsv`.

## Unchanged

Thresholds, sampling-only results, allele-frequency comparisons with SLiM, SLiM false-positive rates (4.9–5.5% on SLiM, 4.9–5.7% on SPOT), mean power differences with SLiM (−0.9 to +1.0 points), SD of power for the trend test (1.5–2.2 and 0.6–1.1 points), false-positive rate with rescaling (4.7–5.2%), and every trend-test result of the case studies quoted in the text (reference power 53%, s for 80% power, power against Ne, depth, human ancient DNA, pathogen).

## Not rerun

- `theory_spot.js` / `theory_check.py`: the exact power without drift is computed with a strict ">" at the threshold; to compare with SPOT's new rule, the exact calculation would need the same tie weight.
- `examples/guide_examples.js` (numbers in `SIMULATION_GUIDE.md`) and `rescale_limits.js`.
- `index.html` keeps its own copy of the simulation code, with the old rule.
- The `.md` write-ups in `validation/` still quote the old numbers.
