# Possible improvements

Ideas that have been discussed but are not implemented yet. SPOT is meant to stay simple, fast and general; for anything beyond its scope, a flexible forward simulator such as [SLiM](https://messerlab.org/slim/) or simuPOP is the better tool. Items here will be considered depending on user feedback.

## Model

- **More flexible demographic histories.** Ne can currently be constant, change exponentially between the first and last sampling points, or be set separately for each interval between sampling points. Possible additions: change points that do not coincide with sampling dates, gradual (rather than abrupt) recovery after a bottleneck, and random fluctuations in population size.
- **Overlapping generations and age structure**, for long-lived species where individuals of different ages breed together.
- **Selection that changes over time**, for example selection that stops, reverses or fluctuates.
- **Population structure**: several connected populations instead of one population with an optional migrant source.
- **Linked selection** (hitchhiking, background selection) and more than one locus.
- **Replicate populations (1 to 100)**: several independent populations starting from the same allele frequency and evolving in parallel (for example several regions or lineages sampled through time), with the evidence combined across replicates, for example by summing their scores. The drift-aware threshold works the same way, and run time grows in proportion to the number of replicates.
- **Sampling a period rather than a date**: spreading the individuals of a time point over a window of years instead of pooling them at one representative year, or giving each individual its own date, read from a file (CSV) rather than typed by hand.

## Sequencing and genotype calls

- **Allele-specific errors**: post-mortem damage (C→T and G→A changes at transitions) and reference bias, which push the frequency in one direction and can create a spurious trend when they differ between dates.
- **Genotype likelihoods**: analysing the reads through genotype likelihoods instead of one called genotype (or one read) per individual.

## Statistics

- Other test statistics (for example likelihood-based tests of selection from time series) next to the trend test.

## Interface

- **Sample-size planner**: columns for the selection coefficient centred on the chosen value (from s/4 to 4s).
- **Power curve**: a logarithmic axis, for example for the selection coefficient.

## Software

- **Continuous integration**: the checks of the `tests` folder (the application gives the same numbers as `sim.js`, also in a real browser) run automatically on GitHub at each change.
- **Browsers**: checking that the same random seed gives the same results in Chrome, Firefox and Safari.
- **Command-line and HPC use**: a small `run.js` wrapper that reads parameters from a JSON file and writes results, so that thousands of settings can be run as an array job on a cluster (for example with SLURM).
- **Python or R package**: a port of the simulation core, if wider use or new features require a full rewrite.
