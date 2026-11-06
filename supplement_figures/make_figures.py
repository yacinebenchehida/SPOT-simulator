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


# =============== Figure S1: agreement with theory (Sections S3.2 and S3.3) ===============
# Selection without drift (Section S3.1) agrees with Equation (7) to rounding error and is not shown.
D = json.load(open('validation/theory_spot.json'))
F = Fig(1060, 420)
# (a) power without drift: SPOT against the exact calculation
TP = rd('validation/theory_power.tsv')
b = Ax(F, 75, 50, 300, 300, (0, 1), (0, 1), title='Power without drift', label='a')
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
# (b) drift: proportion of genetic diversity LOST (as reported by SPOT), all nine scenarios of Table S4
DR = rd('validation/theory_drift.tsv')
c = Ax(F, 500, 50, 330, 300, (0, 90), (0, 0.7), title='Drift: genetic diversity lost', label='b')
c.frame([0, 15, 30, 45, 60, 75, 90], [0, .1, .2, .3, .4, .5, .6, .7], 'generation', 'proportion of genetic diversity lost')
show = ['haploid, N = 100', 'diploid, Nf = 100, Nm = 20', 'diploid, Nf = Nm = 50', 'haplodiploid, Nf = Nm = 100',
        'haplodiploid, Nf = 300, Nm = 50', 'diploid, Ne per interval 1,000 / 50 / 1,000', 'diploid, Nf = Nm = 500, decline to 10%',
        'haploid, N = 1,000', 'diploid, Nf = Nm = 500']
short = {'haploid, N = 100': 'haploid, N = 100', 'diploid, Nf = 100, Nm = 20': 'diploid, Nf = 100, Nm = 20',
         'diploid, Nf = Nm = 50': 'diploid, Nf = Nm = 50', 'haplodiploid, Nf = Nm = 100': 'haplodiploid, Nf = Nm = 100',
         'haplodiploid, Nf = 300, Nm = 50': 'haplodiploid, Nf = 300, Nm = 50',
         'diploid, Ne per interval 1,000 / 50 / 1,000': 'diploid, Ne 1,000 / 50 / 1,000',
         'diploid, Nf = Nm = 500, decline to 10%': 'diploid, decline to 10%', 'haploid, N = 1,000': 'haploid, N = 1,000',
         'diploid, Nf = Nm = 500': 'diploid, Nf = Nm = 500'}
assert set(show) == {r['scenario'] for r in DR}
ends = []
for lab in show:
    rows = [r for r in DR if r['scenario'] == lab]
    assert len(rows) == 91, lab
    sysname = lab.split(',')[0]
    col = SYS[sysname]
    c.line([int(r['generation']) for r in rows], [1 - fl(r['H_ratio_theory']) for r in rows], col, 1.4)
    g = rows[::6]
    c.points([int(r['generation']) for r in g], [1 - fl(r['H_ratio_simulated']) for r in g], col, SH[sysname], 2.8, True)
    ends.append([1 - fl(rows[-1]['H_ratio_theory']), short[lab], col])
# direct labels at the end of each curve, pushed apart vertically when they would overlap
ends.sort(key=lambda e: e[0])
ys = [c.Y(v) for v, _, _ in ends]
for _ in range(50):
    for i in range(1, len(ys)):
        if ys[i - 1] - ys[i] < 12: ys[i] = ys[i - 1] - 12
for (v, txt, col), yy in zip(ends, ys):
    F.text(c.x + c.w + 8, yy + 4, txt, 11, 'start', col)
c.legend([('lines: Equations (5), (6)', COL['muted'], 'line'), ('symbols: SPOT', COL['muted'], 'oh')], c.x + 8, c.y + 16)
print('FigS1', F.save(OUT + 'FigS1_theory.svg'))

# =============== Figure S2: comparison with SLiM ===============
NAME = {'diploid_N1000_s0.02': 'diploid, Ne = 1,000, s = 0.02', 'diploid_N1000_s0.05_rec': 'diploid, Ne = 1,000, s = 0.05, h = 0',
        'diploid_N200_s0.05': 'diploid, Ne = 200, s = 0.05', 'haploid_N1000_s0.02': 'haploid, N = 1,000, s = 0.02',
        'haplodiploid_N1000_s0.02': 'haplodiploid, Ne = 750, s = 0.02', 'haplodiploid_N2000_s0.05': 'haplodiploid, Ne = 1,500, s = 0.05'}
R = rd('validation/slim_power_runs.tsv')
scen = list(dict.fromkeys(r['scenario'] for r in R))
assert set(scen) == set(NAME)
F = Fig(1170, 830)
a = Ax(F, 335, 50, 330, 330, (0.15, 0.8), (0.15, 0.8), title='Power in each of the 20 runs', label='a')   # panel a and its legend centred over panels b-d
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
F = Fig(975, 820)
a = Ax(F, 75, 50, 300, 300, (0, 1), (0, 1), title='Power of the trend test', label='a')
a.frame([0, .2, .4, .6, .8, 1], [0, .2, .4, .6, .8, 1], 'every generation simulated', 'with rescaling', xfmt=pct, yfmt=pct)
xs = [i / 100 for i in range(101)]
se = [2 * math.sqrt(2 * x * (1 - x) / 1000) for x in xs]
a.band(xs, [x - e for x, e in zip(xs, se)], [x + e for x, e in zip(xs, se)], COL['grey'], .25)
a.line([0, 1], [0, 1], COL['muted'], 1, '4 4')
for r in Rs:
    a.points([fl(r['power_exact'])], [fl(r['power_rescaled'])], SYS[r['system']], SH[r['system']], 3.4)
a.legend([(k, v, SH[k]) for k, v in SYS.items()] + [('dashed: 1:1', COL['muted'], 'dash')], a.x + 10, a.y + 18)
# (b) run time against the length of the series, with and without rescaling
b = Ax(F, 640, 50, 300, 300, (3000, 100000), (0.1, 50), xlog=True, ylog=True, title='Run time', label='b')
b.frame([5000, 20000, 60000], [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50], 'generations between the first and last dates',
        'run time (s)', xfmt=lambda v: f'{v:,.0f}', yfmt=lambda v: f'{v:g}')
med = []
for G in ('5000', '20000', '60000'):
    rows = [r for r in Tm if r['G'] == G]
    for j, r in enumerate(rows):
        jit = 1 + 0.035 * ((j % 7) - 3) / 3
        b.points([fl(G) * 0.93 * jit], [fl(r['ms_exact']) / 1000], COL['ink'], 'o', 2.6, True)
        b.points([fl(G) * 1.07 * jit], [fl(r['ms_rescaled']) / 1000], COL['blue'], 'o', 2.6, True)
    lams = sorted({int(r['lambda']) for r in rows})
    F.text(b.X(fl(G)), b.Y(0.1) - 8, 'λ = ' + ' or '.join(str(l) for l in lams), 11, 'middle', COL['muted'])
    med.append((G, st.median(fl(r['ms_exact']) for r in rows) / 1000, st.median(fl(r['ms_rescaled']) for r in rows) / 1000))
b.legend([('every generation simulated', COL['ink'], 'oh'), ('with rescaling', COL['blue'], 'oh')], b.x + 10, b.y + 18)
# (c) limits of rescaling: lambda forced from 2 to 200 (validation/rescale_limits.tsv)
L = rd('validation/rescale_limits.tsv')
c = Ax(F, 240, 460, 360, 300, (1.5, 280), (-0.12, 0.12), xlog=True, title='Rescaling forced beyond its limits', label='c')
c.frame([2, 5, 10, 20, 50, 100, 200], [-.1, -.05, 0, .05, .1], 'rescaling factor λ (forced)',
        'power with rescaling minus without (points)', xfmt=lambda v: f'{v:g}', yfmt=lambda v: f'{v * 100:+.0f}' if v else '0')
c.line([1.5, 280], [0, 0], COL['muted'], 1, '4 4')
DES = [('haploid N=1e4 s=5e-4 p0=0.1', 'haploid, N = 10,000', 'haploid', 1e4, 5e-4),
       ('haploid N=1e3 s=2e-3 p0=0.1', 'haploid, N = 1,000', 'haploid', 1e3, 2e-3),
       ('diploid N=1e5 s=5e-4 p0=0.01', 'diploid, Nf = Nm = 100,000', 'diploid', 1e5, 5e-4),
       ('diploid N=1e3 s=5e-4 p0=0.1', 'diploid, Nf = Nm = 1,000', 'diploid', 1e3, 5e-4)]
assert {d[0] for d in DES} == {r['design'] for r in L}
CL = [COL['blue'], COL['purple'], COL['green'], COL['amber']]
for k, (key, lab, sysname, N, s0) in enumerate(DES):
    lam_max = min(math.floor(math.log(1.025) / math.log(1 + s0)), math.floor(N / 200))   # SPOT's own limits (sim.js, rescale)
    for r in [x for x in L if x['design'] == key]:
        lam = fl(r['lambda']); d = fl(r['difference']) / 100; e = 2 * fl(r['se_difference']) / 100
        xx = lam * (0.88, 0.96, 1.04, 1.12)[k]
        c.errbar([xx], [d - e], [d + e], CL[k], 1.1)
        c.points([xx], [d], CL[k], SH[sysname], 3.6, lam > lam_max)
c.legend([(lab, CL[k], SH[sysname]) for k, (_, lab, sysname, _, _) in enumerate(DES)] +
         [('filled: within SPOT\'s limits', COL['muted'], 'o'), ('hollow: beyond them', COL['muted'], 'oh'),
          ('bars: ±2 standard errors', COL['muted'], 'line')], c.x + c.w + 18, c.y + 16)
print('FigS4', F.save(OUT + 'FigS4_rescaling.svg'), 'median run times (s):', [(g, round(x, 2), round(y, 2)) for g, x, y in med])

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
