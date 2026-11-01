"""Regression of SLiM results on SPOT results (allele frequencies, power, false-positive rates).
Inputs: validation/slim_results.tsv (slim_compare.js), slim_power_v2.tsv and slim_fpr_v2.tsv (slim_compare_v2.js).
Run from the repository root: python3 validation/slim_regression.py
 -> validation/slim_regression.tsv (slope, intercept, R² with 95% confidence intervals) and validation/slim_regression.svg
Ordinary least squares of the SLiM value on the SPOT value; agreement means slope 1, intercept 0."""
import csv, sys, math
sys.path.insert(0, 'validation')
from svgplot import Fig, Ax, COL

rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))
SH = {'haploid': 'o', 'diploid': 's', 'haplodiploid': '^'}
sysof = lambda name: name.split('_')[0]
T975 = {8: 2.306, 31: 2.040, 4: 2.776}   # Student t quantiles for n - 2 degrees of freedom

def ols(x, y):
    n = len(x); mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x); sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    b = sxy / sxx; a = my - b * mx
    res = [yy - a - b * xx for xx, yy in zip(x, y)]; sse = sum(r * r for r in res); syy = sum((v - my) ** 2 for v in y)
    s2 = sse / (n - 2); t = T975[n - 2]
    seb = math.sqrt(s2 / sxx); sea = math.sqrt(s2 * (1 / n + mx * mx / sxx))
    return dict(n=n, slope=b, slope_lo=b - t * seb, slope_hi=b + t * seb, intercept=a, int_lo=a - t * sea, int_hi=a + t * sea,
                r2=1 - sse / syy, max_abs_residual=max(abs(r) for r in res))

F_ = rd('validation/slim_results.tsv')
fa = [(float(r['mean_spot']), float(r['mean_slim']), sysof(r['scenario']), r['run']) for r in F_]
P_ = rd('validation/slim_power_v2.tsv')
pb = [(float(r['mean_power_spot']), float(r['mean_power_slim']), r['system'], 'selected') for r in P_]
R_ = rd('validation/slim_fpr_v2.tsv')
fc = [(float(r['fpr_spot_threshold_on_spot']), float(r['fpr_spot_threshold_on_slim']), r['system'], 'spot') for r in R_] + \
     [(float(r['fpr_sampling_only_on_spot']), float(r['fpr_sampling_only_on_slim']), r['system'], 'naive') for r in R_]

out = [['panel', 'quantity', 'n', 'slope', 'slope_95lo', 'slope_95hi', 'intercept', 'intercept_95lo', 'intercept_95hi', 'R2', 'max_abs_residual']]
F = Fig(1060, 410)
for k, (data, lim, ticks, title, lab, q) in enumerate([
        (fa, (0, 0.5), [0, .1, .2, .3, .4, .5], 'Mean allele frequency', 'a', 'mean frequency at generations 20, 40, 90 (11 scenarios)'),
        (pb, (0, 1), [0, .2, .4, .6, .8, 1], 'Power', 'b', 'mean power over 20 runs (6 scenarios with selection)'),
        (fc, (0, 0.25), [0, .05, .1, .15, .2, .25], 'False-positive rate', 'c', 'false-positive rate, both thresholds (5 neutral scenarios)')]):
    x, y = [d[0] for d in data], [d[1] for d in data]
    r = ols(x, y)
    out.append([lab, q, r['n']] + [f"{r[c]:.4f}" for c in ('slope', 'slope_lo', 'slope_hi', 'intercept', 'int_lo', 'int_hi', 'r2', 'max_abs_residual')])
    ax = Ax(F, 75 + k * 340, 50, 250, 250, lim, lim, title=title, label=lab)
    fmt = (lambda v: f'{v:.0%}') if lab != 'a' else (lambda v: f'{v:g}')
    ax.frame(ticks, ticks, 'SPOT', 'SLiM', xfmt=fmt, yfmt=fmt)
    ax.line([lim[0], lim[1]], [lim[0], lim[1]], COL['muted'], 1.2, '5 4')
    for xx, yy, s_, g in data:
        col = {'selected': COL['blue'], 'neutral': COL['muted'], 'spot': COL['blue'], 'naive': COL['red']}[g]
        ax.points([xx], [yy], col, SH[s_], 4, g == 'neutral')
    F.text(ax.x + 8, ax.y + 18, f"slope {r['slope']:.3f}, R² {r['r2']:.3f}", 11, 'start', COL['ink'])
F.text(75, 372, 'Dashed line: 1:1. Shapes: circle haploid, square diploid, triangle haplodiploid. (a) filled, with selection; hollow, neutral. '
       '(b) each program\'s own threshold.', 11, color=COL['muted'])
F.text(75, 388, '(c) blue, threshold from 20,000 SPOT neutral simulations; red, sampling-only threshold; each applied to 10,000 neutral simulations of each program.', 11, color=COL['muted'])
print('slim_regression.svg', F.save('validation/slim_regression.svg'))
with open('validation/slim_regression.tsv', 'w') as f: f.write('\n'.join('\t'.join(map(str, r)) for r in out) + '\n')
for r in out: print('\t'.join(map(str, r)))
