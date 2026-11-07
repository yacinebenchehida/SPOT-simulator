# SPOT: Selection Power Over Time

**▶ Run it in the browser: <https://yacinebenchehida.github.io/SPOT-simulator/>**

**⬇ Use it offline:** the application is one self-contained file. Save the online page (File → Save page as), or open [`index.html`](index.html) here and click "Download raw file". Then open the file in any browser; no installation and no internet connection are needed.

SPOT estimates the **power to detect selection at a candidate locus** from allele frequencies sampled at two or more dates: museum specimens and a modern sample, an ancient DNA series, or any time series of population samples. It simulates the population forward in time with selection, genetic drift and, optionally, migration and changes in population size, for haploid, diploid or haplodiploid inheritance. It then samples individuals and passes them through low-coverage sequencing and genotype calling, as in a real study.

Power is measured against a **drift-aware threshold**: the detection threshold is set from simulations without selection that include the drift, sampling and genotyping of the design, not from the standard normal distribution. Everything runs locally in the browser, without installation.

## What power means

Power is the probability that a sampling design detects selection when selection is acting. A power of 80% means that selection of the chosen strength is detected in 8 studies out of 10 and missed in 2. SPOT estimates it in five steps:

1. **Simulate the population many times**: by default 1,000 times with selection and 2,000 times without (neutral simulations).
2. **Sample and genotype** individuals at each date, with the sample sizes, read depths and error rates of the design.
3. **Score** each simulation with Z, the Cochran–Armitage trend score: the larger Z, the clearer the rise of the allele over time.
4. **Set the threshold** so that only 5% of the neutral simulations exceed it (the significance level α, divided by the number of candidate loci L).
5. **Count**: power is the proportion of selection simulations above the threshold.

A standard test ignores drift and sets the threshold too low, so it often finds selection where there is none. SPOT shows this as the false-positive rate of the sampling-only threshold.

## Quick start

1. Open the [online version](https://yacinebenchehida.github.io/SPOT-simulator/), or the downloaded file.
2. Describe the design: inheritance system, years per generation, sampling years and number of individuals at each date, sequencing depth, number of candidate loci.
3. Explore what is unknown (Ne, starting frequency p₀, selection coefficient s, dominance, onset of selection) with the **power curve** and the **sample-size planner**, rather than guessing one value.
4. Read the main tile: power of the trend test with the drift-aware threshold. For a reliable estimate, use 5,000 selection simulations.

## Documentation

- **[SPOT manual](SIMULATION_GUIDE.md)**: how the simulations work, step by step with worked examples, and every setting and output of the application.
- **Article**: Ben Chehida et al., in preparation. Its Supporting Information describes the model in full and all validation analyses.

## Using the simulation code in Node.js

`sim.js` holds the simulation functions of the application, with no interface code; the page runs exactly the same code (checked by `validation/engine_sync_check.js`).

```js
const S = require('./sim.js');
const P = {
  system: 'diploid',                  // 'haploid', 'diploid' or 'haplodiploid'
  genTime: 1,                         // years per generation
  Nf: 750, Nm: 750,                   // breeding females and males per generation
  decline: 1,                         // Ne at the last date relative to the first
  mig: 0, pMig: -1,                   // migrants per generation; their allele frequency (-1 = p0)
  p0: 0.1, s: 0.03, h: 0.5, sMale: 1, // starting frequency, selection, dominance, selection in males (multiple of s)
  onset: 1935,                        // year selection starts
  years: [1935, 1955, 1975, 2025],    // sampling years
  n: [20, 20, 20, 40],                // individuals per date
  depthEarly: 10, depthLast: 15,      // mean depth, earlier dates and last date
  err: 0.005, errLast: 0.005, minDP: 3, // per-read error rates and minimum depth
  kSib: 1, maleFrac: 0,               // individuals per family (1 = unrelated); proportion of males sampled
  perfect: false, noDrift: false,     // genotypes from simulated reads; drift on
  neByBlock: false, neBlocks: null,   // or Ne for each interval between dates
  alpha: 0.05, m: 1,                  // significance level; m = number of candidate loci L
  reps: 1000, seed: 1                 // selection simulations (neutral: at least twice as many); random seed
};
const r = S.runPower(P);
console.log(r.trend.powerEmp, r.two.powerEmp);  // power, trend and first-versus-last tests, drift-aware threshold
console.log(r.trend.crit, r.trend.fprNaive);    // drift-aware threshold; false-positive rate of the sampling-only threshold
console.log(1 - S.hetRetained(P), S.neFor(P));  // genetic diversity expected to be lost; Ne
```

`runPower` returns, for the trend test (`trend`) and the first-versus-last test (`two`): the threshold (`crit`), the power with the drift-aware threshold (`powerEmp`) and with the sampling-only threshold (`powerNaive`), and the false-positive rate of the sampling-only threshold (`fprNaive`). It also returns the genotype-call summaries (`hetAcc`, `homAcc`, `hapAcc`, `missRate`) and the proportion of neutral simulations in which drift lost or fixed the allele by the last date (`neutralFixed`). `hetRetained(params)` gives the expected proportion of neutral genetic diversity kept between the first and the last date, and `neFor(params)` the effective population size.

## Reproducing the article

Every number, table and figure of the article comes from a script in this repository. Each script sets its own random seeds and writes a tab-separated file; figures are drawn from these files.

| Result | Script | Results file |
|---|---|---|
| Figure 1 | `figures/make_figure1.js`, `figures/assemble_figure1.py` | `figures/Figure1.png` |
| Figure 2; Tables S1–S3 | `validation/false_positive_rate.js` (letters A, B, C and G); `validation/fpr_figure.py` | `fpr_A.tsv`, `fpr_B.tsv`, `fpr_C.tsv`, `fpr_C_gaussian.tsv` |
| Figure 3; Table S9 (museum insect case study) | `worked_examples/worked_examples.js` (letters A, B, C), `validation/depth_case.js`; `worked_examples/figures.py` | `bee_s_p0.tsv`, `bee_allocation.tsv`, `bee_ne.tsv`, `depth_case.tsv` |
| Example without drift (Section 3.1) | `validation/nodrift_example.js` | `nodrift_example.tsv` |
| Agreement with theory (S3; Table S4, Figure S1) | `validation/theory_spot.js`, `validation/theory_check.py`, `validation/sibling_path_check.js`, `validation/migration_diversity_check.js` | `theory_spot.json`, `theory_power.tsv`, `theory_drift.tsv`, `sibling_path_check.tsv`, `migration_diversity.tsv` |
| Comparison with SLiM (S4; Tables S5–S6, Figure S2) | `validation/slim_check.slim`, `slim_runs.sh`, `slim_runs_v2.sh`, `slim_compare.js`, `slim_compare_v2.js`, `slim_power_spot_threshold.js`, `slim_regression.py` | `slim_*.tsv` |
| Precision and speed (S5; Table S7, Figure S3) | `validation/monte_carlo_precision.js`, `validation/speed.js` | `mc_precision.tsv`, `speed.tsv` |
| Rescaling (S6; Table S8, Figure S4) | `validation/rescale_validation.js`, `rescale_hd_repeat.js`, `rescale_fpr.js`, `rescale_limits.js` | `rescale_*.tsv` |
| Bacterial pathogen (S6.4; Figure S5) | `worked_examples/worked_examples.js` (letter E) | `pathogen.tsv` |
| Human ancient DNA (S7; Table S10, Figure S6) | `worked_examples/worked_examples.js` (letters F, G) | `human_depth.tsv`, `human_budget.tsv` |
| Supplement figures | `supplement_figures/make_figures.py` | `supplement_figures/` |

## How to cite

Please cite the version of SPOT used (shown at the bottom of the page and written into saved settings and results files):

> Ben Chehida, Y. et al. (2026). SPOT: Selection Power Over Time, version 0.9.2. https://github.com/yacinebenchehida/SPOT-simulator

GitHub's "Cite this repository" button gives the same reference in other formats (`CITATION.cff`). The article reference will be added on publication.

## Licence

MIT; see [`LICENSE`](LICENSE).

## Files

| File or folder | Content |
|---|---|
| `index.html` | The application: one self-contained page. |
| `sim.js` | The simulation code, also usable in Node.js. |
| `SIMULATION_GUIDE.md` | The manual. |
| `guide_scripts/` | Scripts that reproduce the numbers of the manual: R scripts for examples 3, 4 and 6, written step by step for teaching; a Node.js script for the others. |
| `figures/` | Figure 1 of the article and its scripts. |
| `validation/` | Validation analyses (scripts and results files), and Figure 2 of the article. |
| `worked_examples/` | Case studies of the article (scripts and results files), and Figure 3. |
| `supplement_figures/` | Figures of the Supporting Information and their script. |
| `LICENSE`, `CITATION.cff` | Licence and citation information. |