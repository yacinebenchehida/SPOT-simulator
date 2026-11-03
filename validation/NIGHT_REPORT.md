# Overnight report (branch validation-v2)

All six tasks were done. Nothing was skipped. `main`, `index.html` and the live app were not touched.

## Commits (all pushed to validation-v2)

| Commit | Task |
|---|---|
| `7d4faa2` | 1. Theory comparison with the tie rule |
| `d354611` | 2. Sibling path check (`P.famPath` in `sim.js`, `validation/sibling_path_check.js`) |
| `c9d966e` | 3 and 4. Figure 3: colours by p₀; panels b and c with three p₀ (one commit: same files) |
| `2cf20e4` | 5. Supporting Information figures (`supplement_figures/`) |
| (this file) | 6. Report |

## 1. Theory comparison with the tie rule

`theory_check.py` now applies the tie rule of `sim.js`: values of Z equal to the threshold count with weight γ (relative tolerance 1e-9, as `TIE_TOL`).
- At the exact 95% quantile, γ is exact, so the exact size is exactly 5%.
- At SPOT's own threshold, SPOT's γ is used. `theory_spot.js` now stores γ and the number of ties.

Outputs regenerated: `theory_power.tsv`, `theory_check.svg`. The legend now reads "neutral-simulation threshold" and "sampling-only threshold".

Power without drift (perfect genotypes), old → new (`theory_check.py`):

| Quantity | Old | New |
|---|---|---|
| Exact size of SPOT's thresholds: mean (range) | 4.39% (1.75–5.61%) | 5.01% (4.43–5.73%) |
| Exact size at the exact 95% quantile: mean (range) | 4.45% (3.21–5.00%) | 5.00% (5.00–5.00%) |
| SPOT power − exact power at SPOT's threshold (neutral-simulation threshold, 120 comparisons): mean; max abs. diff.; abs(z) > 2; RMS z | −0.00 pt; 1.4 pt; 2; 0.89 | +0.00 pt; 1.4 pt; 2; 0.89 |
| SPOT power − exact power at the exact quantile: mean; max abs. diff. | −0.05 pt; 5.9 pt | −0.04 pt; 4.7 pt |
| Sampling-only threshold (120 comparisons) | unchanged | unchanged |

- SPOT's own power without drift changed in 63 of 240 rows of `theory_power.tsv`, by at most 5.9 points.
- The drift (heterozygosity) and selection parts are identical.
- `theory_spot.json` and `theory_drift.tsv` are written by the scripts but had never been committed; they are now in the repository, so that every number of the theory check can be traced.
- `validation/theory.md` still quotes the old numbers.

## 2. Sibling path check

- **Flag:** `P.famPath` in `sim.js` sends k = 1 through the family branch: a mother and a father are drawn from the offspring after selection, then one offspring per family. The default is unchanged: the reference case still gives 53.00% (trend) and 50.08% (first-versus-last).
- **Design:** `validation/sibling_path_check.js` compares the two procedures in diploids and haplodiploids, at Ne = 300, 1,500 and 100,000, with s = 0.02 and 0.05. Reference design, p₀ = 0.1, perfect genotypes, 20,000 selection and 40,000 neutral simulations, 3 seeds; each run sets its own threshold.
- **Output:** `validation/sibling_path_check.tsv`; the run takes about 80 s.

Results:

| Test | Differences (family − default) | Mean | abs(z) > 2 |
|---|---|---|---|
| Trend | −1.4 to +1.0 points | +0.03 | 2 of 33 (max 2.83) |
| First-versus-last | −1.6 to +1.8 points | +0.02 | 4 of 33 (max 3.52) |

Three rows per test have power 100% under both procedures.

The SE in the table is binomial only. It ignores the variation of the threshold between runs, so a few abs(z) > 2 are expected. There is no systematic difference: the two procedures give the same power within Monte Carlo error. The earlier quick test gave −1.1 to +1.7 points.

## 3. Figure 3, colours

`worked_examples/figures.py` now uses the same colour for each p₀ in panels a, b, c and d: 0.01 blue, 0.05 green, 0.1 amber, 0.3 purple. Before, panel d assigned colours by position, so 0.1 was green and 0.3 amber.

Outputs: `bumblebees.svg` and `bumblebees.png` (300 dpi, 3906 × 3906 px).

## 4. Figure 3, panels b and c with three p₀

`worked_examples.js` part C now covers three (p₀, s) pairs; `bee_ne.tsv` gains columns `p0` and `s`. The p₀ = 0.1 rows are identical to the old file.

| p₀ | s | Why this s |
|---|---|---|
| 0.01 | 0.05 | 65% power at Ne = 1,500 in `bee_s_p0.tsv` (67% in `bee_allocation.tsv`) |
| 0.1 | 0.05 | unchanged, as asked |
| 0.3 | 0.015 | see choice below |

Panels b and c: one colour per p₀; the sample-size factors 1, 2 and 5 are shown by line style (solid, dashed, dotted) and marker (circle, square, triangle). The y-axis of c now goes up to 40%.

Numbers for the Results (`worked_examples.js C`; trend test; power at Ne 150 / 300 / 750 / 1,500 / 150,000):

| p₀ (s) | Power ×1 | Power at Ne = 150 with ×2 / ×5 | Sampling-only FPR at Ne = 150, ×1 / ×5 | Max sampling-only FPR |
|---|---|---|---|---|
| 0.01 (0.05) | 20.1 / 30.3 / 50.6 / 65.2 / 98.2% | 17.4 / 20.4% | 3.1 / 4.7% | 11.4% |
| 0.1 (0.05), unchanged | 63.8 / 82.7 / 98.5 / 100 / 100% | 62.6 / 66.7% | 18.7 / 23.8% | 29.0% |
| 0.3 (0.015) | 17.6 / 26.9 / 44.8 / 55.0 / 82.4% | 21.3 / 18.5% | 29.1 / 34.9% | 35.6% |

The conclusion holds for all three p₀: more individuals barely raise power when drift is strong (Ne = 150).

The sampling-only false-positive rate depends strongly on p₀:
- For a rare allele (p₀ = 0.01), it stays close to 5% with the reference samples (1.8–5.8% at all Ne), but rises with larger samples: up to 8.5% with ×2 and 11.4% with ×5, with the maximum at intermediate Ne (750–1,500).
- For p₀ = 0.3, it reaches 35.6%.

## 5. Supporting Information figures

Files are in `supplement_figures/`, as SVG plus PNG at 300 dpi. They are drawn by `make_figures.py` (style of `svgplot.py`), and PNGs are made with `node validation/svg_to_png.js supplement_figures/*.svg`, a new script that renders through the Chromium of Playwright.

| File | Content |
|---|---|
| `FigS2_slim_power_runs` | SPOT against SLiM power in each of the 20 runs, 6 scenarios; filled points = mean |
| `FigS2b_slim_regression` | copy of `validation/slim_regression.svg` |
| `FigS3_precision` | SD of power between 50 runs against R, (a) trend and (b) first-versus-last, with √(P(1 − P)/R) dashed |
| `FigS4_rescaling` | (a) rescaled against full power, 120 cases, by inheritance system, 1:1 line and ±2 SE band; (b) speed-up against λ |
| `FigS5_human_adna` | copy of `worked_examples/human_adna.svg` |
| `FigS6_pathogen` | power against s for the bacterial pathogen (trend, first-versus-last, sampling-only) |

**Speed-up from rescaling:** median 4.5, maximum 33 (manuscript: about 5, up to 30). The run times in `rescale_results.tsv` come from the rerun of 1 November, when four analyses ran in parallel, so the speed-up figures are approximate.

## Choices I made alone

1. **s = 0.015 for p₀ = 0.3** (task 4). No s in the grid of `bee_s_p0.tsv` gives 50–70% power: s = 0.01 gives 33% and s = 0.02 gives 76%. I took s = 0.015, the value already used for p₀ = 0.3 at Ne = 1,500 in the allocation analysis (panel d; 55% power). Panels b, c and d therefore use the same (p₀, s) for p₀ = 0.01 and 0.3. To match d for p₀ = 0.1 as well, s would be 0.02 instead of 0.05, but I kept 0.05 as instructed.
2. **p₀ = 0.1 with s = 0.05 is not intermediate:** its power is about 100% at Ne = 1,500 in panel b. I kept it as asked; you may prefer s = 0.02 (55%), as in panel d.
3. **Tasks 3 and 4 in one commit,** because they change the same files.
4. **Unnumbered regression figure:** I named it `FigS2b`, as part of the SLiM comparison (S2).
5. **New helper `validation/svg_to_png.js`:** no SVG converter was installed.
6. **SE in the sibling check:** binomial only (see section 2).

## Not done or still open

- **Panel 3e** still uses blue and green for the two Ne. These are the colours of p₀ = 0.01 and 0.05 in the other panels, which may confuse readers; consider grey and black, or other colours.
- **Old labels:** `contrasts.svg`, `monte_carlo_precision.svg` and `slim_check.svg` (not in `supplement_figures/`) still say "drift-aware" and "naive".
- **Old numbers in the write-ups:** `validation/theory.md` and the other `.md` files still quote old numbers (see also `ties_summary.tsv`).
- **Manuscript text:** the Results text and the caption of Figure 3 (b, c: "s = 0.05, p₀ = 0.1") must be updated for the three p₀.
