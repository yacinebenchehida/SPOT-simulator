# Figure fixes, 4 November 2026 (branch validation-v2)

Changes made following Yacine's decisions of 4 November (`HANDOFF_figures.md`). `main`, `index.html` and the live app were not touched.

## Commits

| Commit | Task |
|---|---|
| `4f3a40f` | 1, 2. Figure 3: p₀ = 0.1 with s = 0.02 in panels b and c; panel e in dark grey and black |
| `17eaf7e` | 3. Rescaling run times measured alone (`rescale_times.tsv`); Figure S4b uses them |
| `cf1d974` | 4. Notation L; labels "neutral-simulation threshold" and "sampling-only threshold" |
| (this file) | 5. Report |

## 1. Figure 3, panels b and c: p₀ = 0.1 with s = 0.02

The three (p₀, s) pairs are now (0.01, 0.05), (0.1, 0.02) and (0.3, 0.015), the same as in panel d (`worked_examples.js C`, `bee_ne.tsv`). Trend test; ×1, ×2, ×5 = sample sizes 20/20/20/40 multiplied by 1, 2 and 5.

| p₀ (s) | Power ×1 at Ne 150 / 300 / 750 / 1,500 / 150,000 | Power at Ne = 150, ×2 / ×5 | Sampling-only FPR at Ne = 150, ×1 / ×5 | Max sampling-only FPR |
|---|---|---|---|---|
| 0.01 (0.05) | 20.1 / 30.3 / 50.6 / 65.2 / 98.2% | 17.4 / 20.4% | 3.1 / 4.7% | 11.4% |
| 0.1 (0.02), new | 21.6 / 24.8 / 39.4 / 57.9 / 85.0% | 20.0 / 21.0% | 18.5 / 24.9% | 27.1% |
| 0.1 (0.05), old | 63.8 / 82.7 / 98.5 / 100 / 100% | 62.6 / 66.7% | 18.7 / 23.8% | 29.0% |
| 0.3 (0.015) | 17.6 / 26.9 / 44.8 / 55.0 / 82.4% | 21.3 / 18.5% | 29.1 / 34.9% | 35.6% |

- The rows for p₀ = 0.01 and 0.3 are unchanged since the night report.
- For p₀ = 0.1, the sampling-only false-positive rate changes slightly although it does not depend on s. The neutral simulations come from a different part of the random sequence when s changes, so this is Monte Carlo noise only.
- With s = 0.02, all three p₀ have intermediate power at Ne = 1,500 (55–65%). With strong drift (Ne = 150), larger samples barely change power for any p₀ (×5: 20–21%).
- The Results text quoted the old p₀ = 0.1 numbers (64%, 83%, 98.5%; 67% with ×5; sampling-only 19% and 24%); it needs updating.

## 2. Figure 3, panel e

The two effective sizes are now dark grey (Ne = 1,500) and black (Ne = 100,000) instead of blue and green, which are p₀ colours in the other panels. Marker shapes for the error rates are unchanged.

## 3. Rescaling run times, measured alone

`rescale_validation.js times [G]` reruns the two `runPower()` calls of each test case (full, then rescaled) and records only their run times. It was run in three batches (5,000, 20,000 and 60,000 generations) with nothing else running on the machine; the output is `validation/rescale_times.tsv` (120 test cases).

| Series | Cases | λ | Speed-up measured alone: median / max (min) | Old (parallel runs): median / max |
|---|---|---|---|---|
| All | 120 | 3–30 | 4.4 / 28.1 (2.3) | 4.5 / 33.4 |
| 5,000 generations | 54 | 3 | 2.7 / 2.9 | 2.7 / 3.7 |
| 20,000 generations | 54 | 5, 10 | 7.9 / 9.3 | 7.6 / 10.8 |
| 60,000 generations | 12 | 30 | 26.0 / 28.1 | 26.7 / 33.4 |

The manuscript says "median factor of about five, and up to 30". With these times, it is a median of 4.4 and up to 28. Figure S4b now uses these times.

## 4. Notation and labels

- **Figure 2c:** the x-axis reads "per-locus α′ = 0.05 / L". Figure 2 is now also a PNG at 300 dpi (`validation/false_positive_rate.png`).
- **"drift-aware" and "naive" replaced** by "neutral-simulation threshold" and "sampling-only threshold" in:
  - `contrasts.svg` (`figures.py`);
  - `monte_carlo_precision.svg` (`mc_figure.py`);
  - `slim_check.svg` (`slim_figure.py`). The panel c title became "Power (trend test)", because the longer title overlapped panel d; the threshold is named in the footnote instead.
- **No old label is left** in `validation/`, `worked_examples/` or `supplement_figures/`.
- **No other figure** uses the migration rate or the sex subscripts, so nothing else needed changing.

## Files changed

- **Scripts:** `worked_examples/worked_examples.js`, `worked_examples/figures.py`, `validation/rescale_validation.js`, `supplement_figures/make_figures.py`, `validation/fpr_figure.py`, `validation/mc_figure.py`, `validation/slim_figure.py`.
- **Data:** `worked_examples/bee_ne.tsv`, `validation/rescale_times.tsv` (new).
- **Figures:**
  - `worked_examples/bumblebees.svg` / `.png`;
  - `worked_examples/contrasts.svg`;
  - `validation/false_positive_rate.svg` / `.png`;
  - `validation/monte_carlo_precision.svg`;
  - `validation/slim_check.svg`;
  - `supplement_figures/FigS4_rescaling.svg` / `.png`.
