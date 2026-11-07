# SPOT manual

SPOT (Selection Power Over Time) estimates the power to detect selection at a candidate locus from allele frequencies sampled at two or more dates. This manual explains how to use the application and how to read its results. An [explained guide with worked examples](miscellaneous/simulations_explained.md), written for teaching, describes the simulations step by step.

## Contents

1. [Introduction](#1-introduction)
2. [Getting started](#2-getting-started)
3. [Describing a study: the settings](#3-describing-a-study-the-settings)
4. [Using SPOT step by step](#4-using-spot-step-by-step)
5. [Reading and interpreting the results](#5-reading-and-interpreting-the-results)
6. [Questions and answers](#6-questions-and-answers)
7. [Running SPOT from a script](#7-running-spot-from-a-script)

---

## 1. Introduction

### 1.1 What SPOT does

SPOT is designed for planning studies that compare allele frequencies over time: museum specimens and a modern sample, an ancient DNA series, or any time series of population samples. For a candidate locus and a sampling design, it estimates the probability of detecting selection of a given strength (the power), taking into account genetic drift, the number of individuals sampled at each date, sequencing depth and genotyping errors. It can be used to compare designs before sampling or sequencing, for example to decide how many individuals to sequence at each date, or whether a study is worth doing.

### 1.2 How to cite

Please cite SPOT as:

> Ben Chehida, Y. et al. (2026). SPOT: Selection Power Over Time, version 0.9.2. https://github.com/yacinebenchehida/SPOT-simulator

### 1.3 The model in brief

For each run, SPOT:

1. simulates the population forward in time, generation by generation, from the first to the last sampling date, with selection on the favoured allele, genetic drift (Wright–Fisher model) and, if set, migration and changes in population size;
2. samples the chosen number of individuals at each date;
3. simulates sequencing reads for each individual and calls its genotype, with the chosen depth and error rates;
4. computes the score Z of the Cochran–Armitage trend test from the called genotypes;
5. repeats steps 1–4 many times with selection (selection simulations) and without selection (neutral simulations);
6. sets the detection threshold so that only a proportion α′ of the neutral simulations exceed it (the drift-aware threshold), and reports power as the proportion of selection simulations above this threshold.

The equations are given, with their references, in "Model details" at the bottom of the application, and in full in the Supporting Information of the article.

### 1.4 Assumptions and limitations

- One biallelic locus. There is no linkage, so hitchhiking and background selection are not modelled.
- Discrete, non-overlapping generations. For long-lived species with overlapping generations, the generation time and Ne absorb most of the effect.
- Selection is constant in strength from the year it starts and acts on survival before breeding.
- One population, apart from the optional source of migrants: no population structure.
- Population size is constant, changes exponentially between the first and last points, or is set for each interval between sampling points. There is no random fluctuation.
- Families have one mother and one father, so there is no multiple paternity within a family.
- Genotype calling uses a simple per-read error model with a flat prior, with one error rate for the earlier points and one for the last point. Reference bias, mapping quality and damage patterns that depend on the age of each sample are not modelled.
- Periods longer than 2,000 generations are simulated with rescaling: accurate within the tested limits, approximate beyond them (flagged on the page).
- Power is computed for one locus at a time. Testing several candidate loci is handled with the Bonferroni correction (L); genome-wide scans need dedicated methods.

---

## 2. Getting started

SPOT runs in a web browser and needs no installation.

- **Online:** open <https://yacinebenchehida.github.io/SPOT-simulator/>.
- **Offline:** the application is a single file. Save the online page (File → Save page as), or download [`index.html`](index.html) from this repository ("Download raw file"), and open it in any browser. Everything runs on your own computer, and nothing is sent anywhere.

The page has three parts: the **settings** on the left, the **results** on the right (tiles and charts, updated after each change), and, at the bottom, how power is computed, the effect of each setting on power, the model details and references.

The page opens with a diploid example: four sampling dates (1935, 1955, 1975 and 2025, with 20, 20, 20 and 40 individuals), one generation per year, Ne = 2,000, p₀ = 0.1 and s = 0.03. The "Example scenarios" buttons load other designs (Section 3.6).

---

## 3. Describing a study: the settings

This section lists every setting, with its label on the page, its name in the code (`sim.js`, for scripts; Section 7) and the value the page opens with.

### 3.1 Inheritance and population

| Label on the page | Name in the code | Opens with | Meaning |
|---|---|---|---|
| Haplodiploid / Diploid / Haploid | `system` | Diploid | Inheritance system. **Haplodiploid**: females are diploid and males are haploid, developing from unfertilised eggs, so a son inherits only his mother's allele. **Diploid**: both sexes are diploid. **Haploid**: every individual carries one copy of the locus; there are no sexes, no dominance and no families, so the settings for males, dominance and families are hidden. |
| Years per generation | `genTime` | 1 | Average age of parents when their offspring are born, from 0.0001 (about one hour) to 100 years. The number of generations between the first and last points is (last year − first year) / this value, rounded. Periods longer than 2,000 generations are simulated with rescaling (Section 6, "What happens with very long series?"). |
| Breeding females per generation (Nf) | `Nf` | 1,000 | Number of females that actually reproduce each generation (haploid system: "Breeding individuals per generation (N)"). This is not a census count of all individuals, and it is not Ne itself (see below). In social insects, it is roughly the number of successful reproductive females, not the number of workers. |
| Breeding males per generation (Nm) | `Nm` | 1,000 | Number of males that actually reproduce each generation. Raise it relative to Nf to approximate multiple mating by females. Not used in the haploid system. |
| No drift (infinite population) | `noDrift` | off | Switches drift off. Only useful to reproduce the usual calculation of power without drift (the "Analytic baseline" scenario, Section 3.6). Leave it off for realistic power. |
| Set Ne separately for each interval | `neByBlock`, `neBlocks` | off | When ticked, a column "Ne to next" appears in the sampling design. Each value is Ne from that sampling point to the next, which gives simple histories such as a bottleneck between two dates followed by recovery. Breeding females and males are scaled to give exactly that Ne, keeping their ratio. This replaces the breeding numbers and "Ne at the last point relative to the first". |
| Ne at the last point relative to the first | `decline` | 1 | Exponential change of the breeding numbers, and therefore of Ne, between the first and last sampling years: 1 means a constant size, 0.1 a ten-fold decline, 2 a doubling. |
| Migrants per generation, m | `mig` | 0 | Each generation, this fraction of the breeding adults is replaced by migrants from a source population. Keep it at 0 unless immigration from another population is to be modelled. |
| Allele frequency in the migrant source | `pMig` | Same as p₀ (`-1` in the code) | Allele frequency among migrants. "Same as p₀" is a source where the allele is not under selection. Only used when m > 0. |

**How breeding numbers relate to Ne.** The application converts the breeding numbers into an effective population size, shown under the sliders and in the Ne tile:

- haplodiploid: **Ne = 9·Nf·Nm / (4·Nm + 2·Nf)**; with equal numbers N of each sex, Ne = 1.5 N;
- diploid: **Ne = 4·Nf·Nm / (Nf + Nm)**; with equal numbers N of each sex, Ne = 2 N;
- haploid: **Ne = N**.

For most purposes, Ne is the only thing that matters. If the numbers of each sex are unknown, keep Nf = Nm and treat the pair as a single dial for Ne.

### 3.2 Selection

| Label on the page | Name in the code | Opens with | Meaning |
|---|---|---|---|
| Allele frequency at the first point (p₀) | `p0` | 0.1 | Frequency of the favoured allele at the first sampling year, from 0.001 to 0.5. Values of 0.01 or less represent a rare variant; values around 0.1 represent common standing variation. |
| Selection coefficient s | `s` | 0.03 | Relative fitness advantage per generation, from 0.0001 to 0.6. In diploids, the genotype fitnesses are 1 + s (AA), 1 + hs (Aa) and 1 (aa). Values of 0.01–0.05 are common for natural selection; strong selection, such as insecticide resistance, can be ten times higher. Organisms with very short generations (for example bacteria sampled over centuries) need much smaller values per generation (0.0001–0.001). |
| Dominance in diploid individuals | `h` | Additive (h = 0.5) | Effect of the allele in heterozygotes: 1 = dominant, 0.5 = additive, 0 = recessive. Haploid individuals have no heterozygotes. If unknown, keep 0.5 and check 0 and 1 with the parameter power curve. |
| Selection in males (× s) | `sMale` | 1 | Selection in males as a multiple of s (0 to 1.5). Keep it at 1 when the trait is expressed in both sexes; set 0 if it is not expressed in males. |
| Selection starts in year | `onset` | first sampling year | Calendar year when selection begins. Before it, the allele only drifts. Selection acts from the generation after this year. |

### 3.3 Sampling design

| Label on the page | Name in the code | Opens with | Meaning |
|---|---|---|---|
| Year (sampling design) | `years` | 1935, 1955, 1975, 2025 | Calendar years of the sampling points: 2 to 20 points, added with "+ Add time point" and removed with ×. Years must increase and can be negative (before year 0). Samples collected over a period (for example a decade) can be pooled and given one representative year, such as the median. For long-lived species, an individual's alleles reflect the population at its birth, so the birth year (or an estimate of it) is the better date to enter. |
| Individuals (sampling design) | `n` | 20, 20, 20, 40 | Number of individuals sampled at each point (up to 1,000). The first and last points need at least 1. Intermediate points can be set to 0 to drop them; with only the first and last points sampled, the trend test equals the first-versus-last test. |
| Males in the sample | `maleFrac` | 0 | Probability that each sampled individual is male. In haplodiploids, males are haploid and carry one allele instead of two. In diploids, males carry two alleles and differ from females only if "Selection in males" is not 1. 0 means females only. |
| Sampled individuals per family | `kSib` | 1 | 1 means that all sampled individuals are unrelated. Larger values sample families of this size: brothers and sisters with the same mother and father. Siblings carry correlated alleles (in haplodiploids, sisters share their father's allele exactly), which lowers the information in the sample. Raise it only if related individuals are expected in the sample. |

### 3.4 Sequencing and genotype calls

| Label on the page | Name in the code | Opens with | Meaning |
|---|---|---|---|
| Perfect genotypes (ignore sequencing) | `perfect` | off | True genotypes are observed directly, without reads or errors. |
| Mean depth, earlier points | `depthEarly` | 10 | Mean number of reads at the locus for individuals sampled at every point except the last (1 to 30). The number of reads of each individual is drawn from a Poisson distribution with this mean. |
| Mean depth, last point | `depthLast` | 15 | The same for the last point, so that older material can have a lower coverage. |
| Per-read error rate, earlier points | `err` | 0.5% | Probability that a read shows the wrong allele (sequencing error and damage), for every point except the last. |
| Per-read error rate, last point | `errLast` | 0.5% | The same for the last point. Historical samples often carry more errors, so the earlier points can be given a higher rate. All earlier points share the same rate, whatever their age. |
| Minimum depth to call a genotype | `minDP` | 3 | Individuals with fewer reads at the locus are left out of the analysis, so at a low mean depth fewer individuals enter the test. |

Diploid genotypes are called by maximum likelihood among the three genotypes, with a flat prior. At low depth, heterozygotes are often called homozygous because, by chance, only one allele was read (allelic dropout). Haploid individuals are called by the majority of reads, with ties broken at random.

### 3.5 Test and simulation

| Label on the page | Name in the code | Opens with | Meaning |
|---|---|---|---|
| Significance level α | `alpha` | 0.05 | Significance level of the study (0.05, 0.01 or 0.001). |
| Candidate loci tested, L (each tested at α′ = α/L) | `m` | 1 | Number of candidate loci to be tested. Each is tested at α′ = α/L (Bonferroni correction). |
| Selection simulations (at least twice as many neutral) | `reps` | 1,000 | Number of simulations with selection (500, 1,000, 2,000 or 5,000). At least twice as many neutral simulations are run. More simulations give a more precise estimate, not a different answer: with 1,000, power varies between runs by about 1.5–2 percentage points (one standard deviation; [`validation/monte_carlo_precision.md`](validation/monte_carlo_precision.md)). For a final estimate, we recommend 5,000. |
| Random seed | `seed` | 1 | The same seed and the same settings always give exactly the same results, on any computer. |

**Many candidate loci.** The drift-aware threshold is the value of Z exceeded by a proportion α′ of the neutral simulations, which needs about 10 neutral simulations beyond it. When α′ is small, SPOT therefore runs more neutral simulations: 10/α′, up to 200,000 (for example 20,000 for L = 100 and 200,000 for L = 1,000 at α = 0.05). Runs with more than 20,000 neutral simulations are computed in the background, so the page stays responsive. Only beyond 200,000 (more than 1,000 loci at α = 0.05) is the threshold taken from a normal distribution fitted to the neutral values of Z; the status line then says so. To stay fast, the parameter power curve, the planner and the sampling comparison use at most 20,000 neutral simulations per setting; hollow markers in the power curve show points where the normal fit was needed.

### 3.6 Example scenarios

One-click starting points that show the range of designs SPOT covers. They are illustrative values, not estimates for real species, and every value can be changed afterwards.

| Scenario | Settings |
|---|---|
| **Bumblebees** | Haplodiploid; sampling in 1935, 1955, 1975 and 2025 (20, 20, 20 and 40 females); one generation a year; 1,000 breeding females and 1,000 males; p₀ = 0.01, s = 0.05. |
| **Long-lived mammal** | Diploid; 20 years per generation, so only 7 generations from 1880 to 2020; Ne 5,000, then a bottleneck to 200 after 1950 ("Set Ne separately for each interval"); p₀ = 0.1, s = 0.2. Few generations limit power even under strong selection. |
| **Human ancient DNA** | Diploid; 28 years per generation; six sampling dates from 7000 BCE (entered as −7000) to 2000; 50 ancient individuals per date at 2× and 100 present-day individuals at 20×; per-read error rate 1%, minimum depth 1; Ne 10,000; p₀ = 0.1, s = 0.005. |
| **Bacterial pathogen** | Haploid; genomes sampled from 1850 to 2020; 0.01 years (about 4 days) per generation, so about 17,000 generations, simulated with rescaling (Section 6, "What happens with very long series?"); Ne 100,000; p₀ = 0.05, s = 0.0002. |
| **Large diploid population** | Diploid; one generation a year; 1 million breeding females and males, so drift is weak and the sampling-only threshold is almost right; p₀ = 0.05, s = 0.02. |
| **Analytic baseline** | Diploid, with no drift, perfect genotypes and a dominant allele; p₀ = 0.01, s = 0.05. Only the random choice of individuals remains, which reproduces the usual calculation of power without drift (see "How power is computed" in the application). |

---

## 4. Using SPOT step by step

The numbers below come from the settings the page opens with (Section 2). Because each run uses new random simulations, the application may show values that differ by 1 or 2 percentage points.

### 4.1 Power of a sampling design

1. Enter the design: inheritance system, years per generation, sampling years and number of individuals at each date, sequencing depth and error rates, and the number of candidate loci.
2. Enter plausible values for Ne, p₀ and s (Section 4.5 if they are unknown).
3. Read the main tile.

With the settings the page opens with, power is 69% for the trend test and 67% for the first-versus-last test. The drift-aware threshold is Z = 2.11. With the threshold of a standard test that ignores drift (Z = 1.645), power would appear to be 78%, but that test would call 10% of neutral loci "selected" instead of 5%.

### 4.2 Exploring a parameter: the parameter power curve

In "Parameter power curve", choose the parameter to vary and click "Draw curve". Every other setting stays fixed, and the whole power calculation is rerun at each value.

For example, varying "Individuals at the last point" (the earlier points keep 20 individuals each) gives a power of about 47% with 10 individuals, 60% with 20, 69% with 40, 75% with 100 and 80% with 200: each additional individual at the last point adds less and less power, because the earlier samples and drift limit the information.

Varying Ne shows how much drift matters: power is about 42% with Ne = 500, 69% with Ne = 2,000, 84% with Ne = 10,000 and 87% with Ne = 50,000. The lower panel shows that the false-positive rate of the sampling-only threshold rises as Ne falls (16% at Ne = 500, 5% at Ne = 50,000).

### 4.3 Choosing a sample size: the sample-size planner

In "Sample-size planner", choose the parameter for the columns (for example the selection coefficient s) and click "Run planner". Each row multiplies the number of individuals at every date by the same factor (×0.5 to ×4), and each cell gives the power of the trend test. Cells reaching 80% are outlined.

With the settings the page opens with, doubling every sample (40, 40, 40 and 80 individuals) raises power from 69% to 80%; halving it lowers power to 50%.

### 4.4 Haplodiploids: sampling males or females

For haplodiploid species, the panel "Sampling haploid males vs diploid females" ("Compare sampling designs") gives the power of four designs: females only, males only with the same number of individuals, males only with the same number of gene copies (twice as many males), and half males. Males carry one gene copy instead of two, but their genotype is called almost without error. At a low mean depth (up to about 6×), males and females give about the same power per individual; at higher depth, females give slightly more.

### 4.5 Choosing values that are unknown

**Set these to match the study design:** inheritance system, years per generation, sampling years and number of individuals at each point, sex of the sampled individuals, sequencing depth and error rates, and number of candidate loci.

**Explore these with the parameter power curve and the planner:** Ne (through the breeding numbers, or for each interval for a known demographic history), p₀, s, h and the year selection starts. They are rarely known, and power depends strongly on them.

**Leave these as the page opens unless there is a specific reason:** "No drift" (off), "Ne at the last point relative to the first" (1), "Set Ne separately for each interval" (off), "Migrants per generation, m" (0), "Selection in males" (1), "Sampled individuals per family" (1), minimum depth, α, number of simulations and seed.

**Anchoring Ne with diversity data.** If the change in genome-wide heterozygosity between the first and last sampling points can be measured, it gives a plausible Ne (without migration):

```
H_G/H_0 = (1 − 1/(2·Ne))^G          Ne = 1 / ( 2 · (1 − (H_G/H_0)^(1/G)) )
```

Diversity lost over 100 generations for a range of Ne:

| Ne | 50 | 100 | 250 | 500 | 1,000 | 2,500 | 5,000 | 10,000 |
|---|---|---|---|---|---|---|---|---|
| Diversity lost | 63% | 39% | 18% | 9.5% | 4.9% | 2.0% | 1.0% | 0.5% |

The tile "Genetic diversity expected to be lost between the first and the last point" shows this value for the current settings, so the breeding numbers can be adjusted until it matches observed data.

### 4.6 Saving, sharing and reproducing results

- **Copy link** copies the page address, which contains all the current settings: the link reopens SPOT with exactly the same design.
- **Save settings** and **Load settings** save the settings to a small file and load them again later.
- **Download results (CSV)** saves the settings and the main results. The CSV buttons of the parameter power curve and the planner save the values behind each chart, and the PNG and SVG buttons save the charts as images.
- With the same settings and the same random seed, SPOT gives exactly the same results on any computer.

---

## 5. Reading and interpreting the results

**Tiles**

| Tile | Meaning |
|---|---|
| **Power · trend test across all sampling points, drift-aware threshold** | Main result: the proportion of selection simulations whose score Z, computed from all sampling points, exceeds the drift-aware threshold. Below it, the same for the first-versus-last test. A note appears when, in more than half of the neutral simulations, drift has lost or fixed the allele by the last point. |
| Power · trend test, sampling-only threshold | Power with the fixed threshold of a standard test that ignores drift (1.645 for α′ = 5%), and the false-positive rate of that threshold: the proportion of neutral simulations it calls "selected". A warning appears when this threshold overstates power by 5 points or more. |
| Genetic diversity expected to be lost between the first and the last point | Average over neutral loci, from Ne (formula in "Model details" in the application). With migration, the value with migration is shown next to the value without it. |
| Effective population size Ne | Computed from the breeding numbers (first → last point if the population changes), or the lowest and highest values when Ne is set for each interval. |
| Heterozygotes called correctly | Proportion of true heterozygotes called heterozygous (haploid system: gene copies called correctly), with the proportion of homozygotes called correctly and of sampled individuals below the minimum depth. |

**Charts**

- **Allele trajectories:** 40 of the selection simulations, 40 of the neutral simulations, and the sampling points. The text under the chart gives the expected frequencies without drift.
- **Scores Z with and without selection:** distributions of Z in the neutral simulations and in the selection simulations, with the drift-aware threshold (solid line) and the sampling-only threshold (dashed line). Power is the proportion of the selection distribution beyond the solid line.
- **Parameter power curve:** every setting stays fixed except the parameter chosen, which is stepped through a range of values; at each value, the whole power calculation is rerun (new simulations with and without selection, new threshold). The upper panel shows power (trend test and first-versus-last test, drift-aware threshold); the lower panel shows the false-positive rate of the sampling-only threshold. The vertical line marks the current value.
- **Sample-size planner:** a grid of power (trend test, drift-aware threshold). Rows scale the whole sampling design (×0.5 to ×4); columns vary one parameter (s, h, Ne or p₀). Cells reaching 80% are outlined. Each cell uses up to 1,000 selection simulations.
- **Sampling haploid males vs diploid females:** power for four sampling designs, with the genotype calls shown on hover. Shown only for haplodiploids (in diploids, both sexes carry two alleles).

**Interpreting power.**

- The main result is the power of the trend test with the drift-aware threshold. Studies are usually planned to reach a power of at least 80%.
- If the first-versus-last test gives almost the same power as the trend test, the intermediate dates add little information for this scenario.
- Power is estimated from simulations. With the default 1,000 selection simulations, it varies by about 1.5–2 percentage points between runs; for a final estimate, we recommend 5,000 selection simulations.
- The sampling-only tile shows what a standard test that ignores drift would conclude. When its power is much higher than the power with the drift-aware threshold, ignoring drift would overstate the power of the design, and its false-positive rate shows how often such a test would call neutral loci "selected".
- Power applies to one candidate locus, for the chosen values of Ne, p₀, s, dominance and the year selection starts. When these are uncertain, the parameter power curve and the planner show how robust the design is (Section 4.5).

**Notes and warnings shown on the page.**

- *Drift switched off:* power then reflects sampling only, and the drift-aware threshold is close to the threshold of a standard test.
- *Allele lost or fixed by drift:* shown when drift has lost or fixed the allele by the last date in more than half of the neutral simulations; little variation is then left to test.
- *Hollow markers in the parameter power curve:* too few neutral simulations for the chosen α′, so the threshold was taken from a normal distribution fitted to the neutral values; these points can be far off.
- *Rescaling:* very long series are simulated with rescaling (Section 6); results outside the tested limits are flagged as approximate.

---

## 6. Questions and answers

**Ne is unknown. What should I enter?**
Explore a range of values with the parameter power curve. If the loss of genetic diversity between the first and last dates has been measured, it gives a plausible Ne (Section 4.5).

**Which selection coefficient should I use?**
s is the fitness advantage per generation. Values of 0.01–0.05 are common for natural selection; strong selection, such as insecticide resistance, can be ten times higher. For organisms with very short generations, values per generation are much smaller (0.0001–0.001). Exploring a range with the planner is usually more informative than a single guess.

**Why does power collapse when there are very many generations between the dates?**
Over many generations, drift loses or fixes most neutral alleles, and selected alleles can be fixed before the first date or between dates; the frequencies then no longer change between samples. SPOT shows a note when this happens in most neutral simulations.

**Why is the false-positive rate of the sampling-only threshold lower at low depth?**
Individuals with fewer reads than the minimum depth are left out, so fewer individuals enter the test. Sampling then dominates the chance changes in frequency, and the standard threshold, which accounts for sampling only, is closer to correct.

**Why does the false-positive rate of the sampling-only threshold not rise at small Ne when the allele is rare?**
With a rare starting allele (for example p₀ = 0.01), drift often loses it entirely; the score Z is then 0 and cannot exceed the threshold.

**How many simulations should I use?**
The default 1,000 selection simulations are enough to explore designs. For a final estimate, we recommend 5,000.

**Several candidate loci will be tested. What changes?**
Set "Candidate loci tested, L". Each locus is then tested at α′ = α/L (Bonferroni correction), which raises the threshold and lowers power.

**What happens with very long series?**
Historical DNA of organisms with short generations (for example bacteria sampled over decades to centuries) can span hundreds of thousands of generations. Simulating each of them would take minutes, so SPOT uses the standard rescaling of forward population-genetic simulators (as in SLiM): when the period spans more than 2,000 generations, each simulated generation stands for λ real ones, with

```
selection          s'   = (1 + s)^λ − 1          (same for selection in males)
breeding numbers   N'   = N / λ                   (also Ne for each interval)
migration          m'   = 1 − (1 − m)^λ
sampling times     t'   = t / λ
```

This keeps Ne·s and the total drift over the period unchanged; sampling, genotype calling and the test are not affected. λ is chosen so that about 2,000 generations are simulated, but it is capped so that the rescaled selection stays at or below 0.025 and the rescaled breeding numbers at or above 200. These limits come from comparisons with simulations of every generation ([`validation/rescaling.md`](validation/rescaling.md); Supporting Information S6 of the article). The status line says when rescaling is used. If the caps would require simulating more than 10,000 generations, λ is raised anyway to keep the page fast, and a note marks the result as approximate.

**Can SPOT analyse my data?**
No. SPOT is for planning: it estimates power for a design before the data exist. To test real data in the same spirit, the counterpart of the neutral simulations is an empirical null distribution: the score Z computed for genome-wide SNPs with a similar starting frequency, in the same samples.

**Why do two runs give slightly different results?**
Each run draws new random simulations. With the same settings and the same random seed, the results are exactly the same.

---

## 7. Running SPOT from a script

The application is all that is needed to use SPOT. For many runs (for example a large grid of designs), the simulations can also be run from a script with [Node.js](https://nodejs.org). `sim.js` holds the simulation functions of the application, with no interface code; the page runs exactly the same code (checked by `validation/engine_sync_check.js`). The settings are those of Section 3, under their names in the code.

```js
const S = require('./sim.js');
const P = {
  system: 'diploid',                  // 'haploid', 'diploid' or 'haplodiploid'
  genTime: 1,                         // years per generation
  Nf: 1000, Nm: 1000,                 // breeding females and males per generation
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
  alpha: 0.05, m: 1,                  // significance level; number of candidate loci L (m in the code; migration is mig)
  reps: 1000, seed: 1                 // selection simulations (neutral: at least twice as many); random seed
};
const r = S.runPower(P);
console.log(r.trend.powerEmp, r.two.powerEmp);  // power, trend and first-versus-last tests, drift-aware threshold
console.log(r.trend.crit, r.trend.fprNaive);    // drift-aware threshold; false-positive rate of the sampling-only threshold
console.log(1 - S.hetRetained(P), S.neFor(P));  // genetic diversity expected to be lost; Ne
```

`runPower` returns, for the trend test (`trend`) and the first-versus-last test (`two`): the threshold (`crit`), the power with the drift-aware threshold (`powerEmp`) and with the sampling-only threshold (`powerNaive`), and the false-positive rate of the sampling-only threshold (`fprNaive`). It also returns the genotype-call summaries (`hetAcc`, `homAcc`, `hapAcc`, `missRate`) and the proportion of neutral simulations in which drift lost or fixed the allele by the last date (`neutralFixed`). `hetRetained(params)` gives the expected proportion of neutral genetic diversity kept between the first and the last date, and `neFor(params)` the effective population size.
