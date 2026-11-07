# Analyses and figures of the article

This folder contains everything needed to reproduce the results of the article describing SPOT (Ben Chehida et al., in preparation) and of its Supporting Information. All scripts use the same simulation code as the application (`../sim.js`). Run them from the root of the repository; each one sets its own random seeds and writes a tab-separated results file, from which the figures are drawn.

| Result | Script | Results file |
|---|---|---|
| Figure 1 | `figures/make_figure1.js`, `figures/assemble_figure1.py` | `figures/Figure1.png` |
| Figure 2; Tables S1–S3 | `validation/false_positive_rate.js` (letters A, B, C and G); `validation/fpr_figure.py` | `validation/fpr_A.tsv`, `fpr_B.tsv`, `fpr_C.tsv`, `fpr_C_gaussian.tsv` |
| Figure 3; Table S9 (museum insect case study) | `worked_examples/worked_examples.js` (letters A, B, C), `validation/depth_case.js`; `worked_examples/figures.py` | `worked_examples/bee_s_p0.tsv`, `bee_allocation.tsv`, `bee_ne.tsv`, `depth_case.tsv` |
| Example without drift (Section 3.1) | `validation/nodrift_example.js` | `validation/nodrift_example.tsv` |
| Agreement with theory (S3; Table S4, Figure S1) | `validation/theory_spot.js`, `validation/theory_check.py`, `validation/sibling_path_check.js`, `validation/migration_diversity_check.js` | `validation/theory_spot.json`, `theory_power.tsv`, `theory_drift.tsv`, `sibling_path_check.tsv`, `migration_diversity.tsv` |
| Comparison with SLiM (S4; Tables S5–S6, Figure S2) | `validation/slim_check.slim`, `slim_runs.sh`, `slim_runs_v2.sh`, `slim_compare.js`, `slim_compare_v2.js`, `slim_power_spot_threshold.js`, `slim_regression.py` | `validation/slim_*.tsv`; SLiM outputs in `validation/slim_out/` |
| Precision and speed (S5; Table S7, Figure S3) | `validation/monte_carlo_precision.js`, `validation/speed.js` | `validation/mc_precision.tsv`, `speed.tsv` |
| Rescaling (S6; Table S8, Figure S4) | `validation/rescale_validation.js`, `rescale_hd_repeat.js`, `rescale_fpr.js`, `rescale_limits.js` | `validation/rescale_*.tsv` |
| Bacterial pathogen (S6.4; Figure S5) | `worked_examples/worked_examples.js` (letter E) | `worked_examples/pathogen.tsv` |
| Human ancient DNA (S7; Table S10, Figure S6) | `worked_examples/worked_examples.js` (letters F, G) | `worked_examples/human_depth.tsv`, `human_budget.tsv` |
| Supplement figures S1–S6 | `supplement_figures/make_figures.py` | `supplement_figures/` |

Each script starts with a short comment giving the exact command. For example:

```
node paper/validation/false_positive_rate.js A > paper/validation/fpr_A.tsv
python3 paper/supplement_figures/make_figures.py
```

The SLiM comparison also needs [SLiM](https://messerlab.org/slim/) (version 5).

`others/` holds files that were used during the development of the article but are no longer needed to reproduce it (see `others/README.md`).
