"""Figures S1-S6 of the Supporting Information, drawn from the result tables in validation/ and worked_examples/.
Run from the repository root:
    python3 supplement_figures/make_figures.py
    node validation/svg_to_png.js supplement_figures/FigS*.svg      (PNG at 300 dpi)
Writes supplement_figures/FigS1_theory.svg ... FigS6_pathogen.svg. As in the main figures, nothing is written under
the panels: what the symbols and lines mean is in the legends, and everything else is in the figure captions."""
import csv, json, math, statistics as st, sys
sys.path.insert(0, 'validation')
from svgplot import Fig, Ax, COL

rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))
fl = float
pct = lambda v: f'{v * 100:.0f}%'
OUT = 'supplement_figures/'
SYS = {'haploid': COL['blue'], 'diploid': COL['green'], 'haplodiploid': COL['amber']}
SH = {'haploid': 'o', 'diploid': 's', 'haplodiploid': '^'}
PAL = [COL['blue'], COL['green'], COL['amber'], COL['purple'], COL['red'], COL['ink']]
TH = {'neutral': COL['blue'], 'sampling': COL['red']}   # the two thresholds, as in Figure 2


def ols(x, y):
    n = len(x); mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x); sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    b = sxy / sxx; a = my - b * mx
    sse = sum((yy - a - b * xx) ** 2 for xx, yy in zip(x, y)); syy = sum((v - my) ** 2 for v in y)
    return b, 1 - sse / syy


# =============== Figure S1: agreement with theory (panels in the order of Section S3) ===============
D = json.load(open('validation/theory_spot.json'))
F = Fig(1200, 420)
# (a) selection without drift, haploids: SPOT against the closed form (Equation 7)
a = Ax(F, 75, 50, 300, 300, (0, 300), (0, 1), title='Selection without drift', label='a')
a.frame([0, 50, 100, 150, 200, 250, 300], [0, .2, .4, .6, .8, 1], 'generation', 'frequency of the favoured allele')
for k, s in enumerate([0.01, 0.05, 0.2]):
    for p0 in (0.01, 0.2):
        r = next(x for x in D['selection'] if x['system'] == 'haploid' and x['s'] == s and x['p0'] == p0)
        th = [p0 * (1 + s) ** t / (1 - p0 + p0 * (1 + s) ** t) for t in range(301)]
        a.line(list(range(0, 301, 2)), th[::2], PAL[k], 1.6)
        g = list(range(0, 301, 10))
        a.points(g, [r['f'][i] for i in g], PAL[k], 'o', 2.6, True)
a.legend([('s = 0.01', PAL[0], 'line'), ('s = 0.05', PAL[1], 'line'), ('s = 0.2', PAL[2], 'line'),
          ('lines: Equation (7)', COL['muted'], 'line'), ('circles: SPOT', COL['muted'], 'oh')], a.X(170), a.Y(0.5))
# (b) power without drift: SPOT against the exact calculation
TP = rd('validation/theory_power.tsv')
b = Ax(F, 490, 50, 300, 300, (0, 1), (0, 1), title='Power without drift', label='b')
b.frame([0, .2, .4, .6, .8, 1], [0, .2, .4, .6, .8, 1], 'exact power', 'power given by SPOT', xfmt=pct, yfmt=pct)
xs = [i / 100 for i in range(101)]
se = [2 * math.sqrt(x * (1 - x) / 5000) for x in xs]
b.band(xs, [x - e for x, e in zip(xs, se)], [x + e for x, e in zip(xs, se)], COL['grey'], .25)
b.line([0, 1], [0, 1], COL['muted'], 1, '4 4')
for r in TP:
    col = TH['sampling'] if r['threshold'].startswith('naive') else TH['neutral']
    first_last = r['statistic'] == 'two' or r['design'] == '2 points'
    b.points([fl(r['power_exact_at_spot_threshold'])], [fl(r['power_spot'])], col, 's' if first_last else 'o', 3.2)
b.legend([('neutral-simulation threshold', TH['neutral'], 'o'), ('sampling-only threshold', TH['sampling'], 'o'),
          ('trend test', COL['muted'], 'o'), ('first-versus-last test', COL['muted'], 's')], b.x + 8, b.y + 16)
# (c) drift: proportion of genetic diversity kept
DR = rd('validation/theory_drift.tsv')
c = Ax(F, 905, 50, 270, 300, (0, 90), (0, 1), title='Drift: genetic diversity kept', label='c')
c.frame([0, 20, 40, 60, 80], [0, .2, .4, .6, .8, 1], 'generation', 'proportion of genetic diversity kept')
show = [('haploid, N = 100', 'haploid, N = 100'), ('diploid, Nf = Nm = 50', 'diploid, Nf = Nm = 50'),
        ('haplodiploid, Nf = Nm = 100', 'haplodiploid, Nf = Nm = 100'),
        ('diploid, Nf = Nm = 500, decline to 10%', 'diploid, decline to 10%'),
        ('diploid, Ne per interval 1,000 / 50 / 1,000', 'diploid, Ne 1,000, 50, 1,000')]
for k, (lab, short) in enumerate(show):
    rows = [r for r in DR if r['scenario'] == lab]
    assert len(rows) == 91, lab
    c.line([int(r['generation']) for r in rows], [fl(r['H_ratio_theory']) for r in rows], PAL[k], 1.6)
    g = rows[::6]
    c.points([int(r['generation']) for r in g], [fl(r['H_ratio_simulated']) for r in g], PAL[k], 'o', 3, True)
c.legend([(short, PAL[k], 'line') for k, (lab, short) in enumerate(show)] +
         [('lines: Equations (5), (6)', COL['muted'], 'line'), ('circles: SPOT', COL['muted'], 'oh')], c.x + 8, c.y + c.h - 112)
print('FigS1', F.save(OUT + 'FigS1_theory.svg'))

# =============== Figure S2: comparison with SLiM ===============
NAME = {'diploid_N1000_s0.02': 'diploid, Ne = 1,000, s = 0.02', 'diploid_N1000_s0.05_rec': 'diploid, Ne = 1,000, s = 0.05, h = 0',
        'diploid_N200_s0.05': 'diploid, Ne = 200, s = 0.05', 'haploid_N1000_s0.02': 'haploid, N = 1,000, s = 0.02',
        'haplodiploid_N1000_s0.02': 'haplodiploid, Ne = 750, s = 0.02', 'haplodiploid_N2000_s0.05': 'haplodiploid, Ne = 1,500, s = 0.05'}
R = rd('validation/slim_power_runs.tsv')
scen = list(dict.fromkeys(r['scenario'] for r in R))
assert set(scen) == set(NAME)
F = Fig(1200, 830)
a = Ax(F, 75, 50, 330, 330, (0.15, 0.8), (0.15, 0.8), title='Power in each of the 20 runs', label='a')
a.frame([.2, .3, .4, .5, .6, .7, .8], [.2, .3, .4, .5, .6, .7, .8], 'power with SLiM', 'power with SPOT', xfmt=pct, yfmt=pct)
a.line([0.15, 0.8], [0.15, 0.8], COL['muted'], 1, '4 4')
for k, s in enumerate(scen):
    r = [x for x in R if x['scenario'] == s]
    a.points([fl(x['power_slim']) for x in r], [fl(x['power_spot']) for x in r], PAL[k], 'o', 3.2, True)
    a.points([st.mean(fl(x['power_slim']) for x in r)], [st.mean(fl(x['power_spot']) for x in r)], PAL[k], 'o', 5.5)
a.legend([(NAME[s], PAL[k], 'o') for k, s in enumerate(scen)] +
         [('hollow: one run', COL['muted'], 'oh'), ('filled: mean of the 20 runs', COL['muted'], 'o'), ('dashed: 1:1', COL['muted'], 'dash')],
         a.x + a.w + 25, a.y + 18)
fa = [(fl(r['mean_spot']), fl(r['mean_slim']), r['scenario'].split('_')[0], r['run']) for r in rd('validation/slim_results.tsv')]
pb = [(fl(r['mean_power_spot']), fl(r['mean_power_slim']), r['system'], 'selected') for r in rd('validation/slim_power_v2.tsv')]
RF = rd('validation/slim_fpr_v2.tsv')
fc = [(fl(r['fpr_spot_threshold_on_spot']), fl(r['fpr_spot_threshold_on_slim']), r['system'], 'neutral-thr') for r in RF] + \
     [(fl(r['fpr_sampling_only_on_spot']), fl(r['fpr_sampling_only_on_slim']), r['system'], 'sampling-thr') for r in RF]
for k, (data, lim, ticks, title, lab, f_) in enumerate([
        (fa, (0, 0.5), [0, .1, .2, .3, .4, .5], 'Mean allele frequency', 'b', lambda v: f'{v:g}'),
        (pb, (0, 1), [0, .2, .4, .6, .8, 1], 'Mean power', 'c', pct),
        (fc, (0, 0.25), [0, .05, .1, .15, .2, .25], 'False-positive rate', 'd', pct)]):
    ax = Ax(F, 75 + k * 390, 490, 280, 280, lim, lim, title=title, label=lab)
    ax.frame(ticks, ticks, 'SPOT', 'SLiM', xfmt=f_, yfmt=f_)
    ax.line(list(lim), list(lim), COL['muted'], 1.2, '5 4')
    for xx, yy, sy, grp in data:
        col = {'selected': COL['ink'], 'neutral': COL['ink'], 'neutral-thr': TH['neutral'], 'sampling-thr': TH['sampling']}[grp]
        ax.points([xx], [yy], col, SH[sy], 4, grp == 'neutral')
    slope, r2 = ols([d[0] for d in data], [d[1] for d in data])
    F.text(ax.x + 8, ax.y + 18, f'slope {slope:.3f}, R² {r2:.3f}', 11, 'start', COL['ink'])
    if lab == 'b':
        ax.legend([('haploid', COL['ink'], 'o'), ('diploid', COL['ink'], 's'), ('haplodiploid', COL['ink'], '^'),
                   ('filled: with selection', COL['ink'], 'o'), ('hollow: without selection', COL['ink'], 'oh')], ax.x + 8, ax.Y(0.415))
    if lab == 'd':
        ax.legend([('neutral-simulation threshold', TH['neutral'], 'o'), ('sampling-only threshold', TH['sampling'], 'o')], ax.X(0.09), ax.Y(0.06))
print('FigS2', F.save(OUT + 'FigS2_slim.svg'))

# =============== Figure S3: variation of power between runs ===============
M = rd('validation/mc_precision.tsv')
sc5 = list(dict.fromkeys(r['scenario'] for r in M))
F = Fig(1200, 420)
for j, (stat, title) in enumerate((('trend', 'Trend test'), ('two', 'First-versus-last test'))):
    ax = Ax(F, 85 + 430 * j, 50, 320, 300, (400, 6500), (0, 0.04), xlog=True, title=title, label='ab'[j])
    ax.frame([500, 1000, 2000, 5000], [0, .01, .02, .03, .04], 'selection simulations R', 'SD of power (percentage points)',
             xfmt=lambda v: f'{v:,.0f}', yfmt=lambda v: f'{v * 100:.0f}')
    for k, s in enumerate(sc5):
        r = sorted([x for x in M if x['scenario'] == s and x['statistic'] == stat], key=lambda x: fl(x['replicates']))
        xs = [fl(x['replicates']) for x in r]
        ax.line(xs, [fl(x['sd_power']) for x in r], PAL[k]); ax.points(xs, [fl(x['sd_power']) for x in r], PAL[k])
        ax.line(xs, [fl(x['sd_binomial_only']) for x in r], PAL[k], 1.2, '4 3')
ax.legend([(f'"{s}"', PAL[k], 'o-') for k, s in enumerate(sc5)] +
          [('solid: observed between 50 runs', COL['muted'], 'line'), ('dashed: binomial √(P(1 − P)/R)', COL['muted'], 'dash')],
          ax.x + ax.w + 25, ax.y + 18)
print('FigS3', F.save(OUT + 'FigS3_precision.svg'))

# =============== Figure S4: rescaling ===============
Rs = rd('validation/rescale_results.tsv')
Tm = rd('validation/rescale_times.tsv')   # run times measured with nothing else running, same order as Rs
assert [(r['G'], r['system'], r['N'], r['s'], r['p0']) for r in Tm] == [(r['G'], r['system'], r['N'], r['s'], r['p0']) for r in Rs]
F = Fig(1000, 420)
a = Ax(F, 75, 50, 330, 330, (0, 1), (0, 1), title='Power of the trend test', label='a')
a.frame([0, .2, .4, .6, .8, 1], [0, .2, .4, .6, .8, 1], 'every generation simulated', 'with rescaling', xfmt=pct, yfmt=pct)
xs = [i / 100 for i in range(101)]
se = [2 * math.sqrt(2 * x * (1 - x) / 1000) for x in xs]
a.band(xs, [x - e for x, e in zip(xs, se)], [x + e for x, e in zip(xs, se)], COL['grey'], .25)
a.line([0, 1], [0, 1], COL['muted'], 1, '4 4')
for r in Rs:
    a.points([fl(r['power_exact'])], [fl(r['power_rescaled'])], SYS[r['system']], SH[r['system']], 3.4)
a.legend([(k, v, SH[k]) for k, v in SYS.items()] + [('dashed: 1:1', COL['muted'], 'dash')], a.x + 10, a.y + 18)
sp = [fl(r['speedup']) for r in Tm]
b = Ax(F, 545, 50, 330, 330, (2, 40), (1, 100), xlog=True, ylog=True, title='Speed-up from rescaling', label='b')
b.frame([2, 5, 10, 20, 40], [1, 2, 5, 10, 20, 50, 100], 'rescaling factor λ', 'run time: every generation / rescaled')
b.line([2, 40], [2, 40], COL['muted'], 1, '4 4')
for r, v in zip(Rs, sp):
    b.points([fl(r['lambda'])], [v], SYS[r['system']], SH[r['system']], 3.4, True)
b.legend([('dashed: speed-up = λ', COL['muted'], 'dash')], b.x + 10, b.y + 18)
print('FigS4', F.save(OUT + 'FigS4_rescaling.svg'), 'median speed-up', round(st.median(sp), 2), 'max', round(max(sp), 1))

# =============== Figure S5: human ancient DNA case study ===============
F = Fig(1000, 420)
HD = rd('worked_examples/human_depth.tsv')
g = Ax(F, 75, 50, 340, 300, (0.0008, 0.13), (0, 1), xlog=True, title='Depth of the ancient samples', label='a')
g.frame([0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1], [0, .2, .4, .6, .8, 1], 'selection coefficient s', 'power (trend test)', yfmt=pct)
g.hline(0.8, COL['muted'])
deps = ['1', '2', '5', '10', 'perfect']
for k, dp in enumerate(deps):
    r = [x for x in HD if x['ancient_depth'] == dp]
    xs = [fl(x['s']) for x in r]
    g.line(xs, [fl(x['power_trend']) for x in r], PAL[k]); g.points(xs, [fl(x['power_trend']) for x in r], PAL[k])
g.legend([((dp + '×') if dp != 'perfect' else 'perfect genotypes', PAL[k], 'o-') for k, dp in enumerate(deps)], g.x + 10, g.y + 18)
HB = rd('worked_examples/human_budget.tsv')
h = Ax(F, 570, 50, 340, 300, (8, 125), (0, 1), xlog=True, title='100 reads per ancient date', label='b')
h.frame([10, 20, 25, 50, 100], [0, .2, .4, .6, .8, 1], 'individuals per ancient date (mean depth = 100 / individuals)', 'power (trend test)', yfmt=pct)
for k, s in enumerate(['0.005', '0.0075']):
    r = [x for x in HB if x['s'] == s]
    xs = [fl(x['individuals_per_ancient_date']) for x in r]
    h.line(xs, [fl(x['power_trend']) for x in r], PAL[k]); h.points(xs, [fl(x['power_trend']) for x in r], PAL[k])
h.legend([(f's = {s}', PAL[k], 'o-') for k, s in enumerate(['0.005', '0.0075'])], h.x + 10, h.y + 18)
print('FigS5', F.save(OUT + 'FigS5_human_adna.svg'))

# =============== Figure S6: bacterial pathogen, power against s (rescaled) ===============
P = rd('worked_examples/pathogen.tsv')
F = Fig(560, 420)
f = Ax(F, 85, 50, 400, 300, (1.5e-5, 7e-4), (0, 1), xlog=True, title='Bacterial pathogen, 17,000 generations', label=None)
f.frame([2e-5, 5e-5, 1e-4, 2e-4, 5e-4], [0, .2, .4, .6, .8, 1], 'selection coefficient s per generation', 'power', xfmt=lambda v: f'{v:.5f}'.rstrip('0'), yfmt=pct)
f.hline(0.8, COL['muted'])
xs = [fl(x['s_per_generation']) for x in P]
f.line(xs, [fl(x['power_trend']) for x in P], COL['blue']); f.points(xs, [fl(x['power_trend']) for x in P], COL['blue'])
f.line(xs, [fl(x['power_first_vs_last']) for x in P], COL['green']); f.points(xs, [fl(x['power_first_vs_last']) for x in P], COL['green'], 's')
f.line(xs, [fl(x['power_naive']) for x in P], COL['red'], 1.2, '4 3')
f.legend([('trend test', COL['blue'], 'o-'), ('first-versus-last test', COL['green'], 's-'),
          ('trend test, sampling-only threshold', COL['red'], 'dash')], f.x + 10, f.y + 18)
print('FigS6', F.save(OUT + 'FigS6_pathogen.svg'))
