"""Figures of the Supporting Information, from result tables in validation/ and worked_examples/.
Run from the repository root:  python3 supplement_figures/make_figures.py  &&  node validation/svg_to_png.js supplement_figures/*.svg
Writes supplement_figures/FigS2_slim_power_runs.svg, FigS3_precision.svg, FigS4_rescaling.svg, FigS6_pathogen.svg;
FigS2b_slim_regression.svg and FigS5_human_adna.svg are copies of validation/slim_regression.svg and
worked_examples/human_adna.svg (made by slim_regression.py and worked_examples/figures.py)."""
import csv, math, shutil, statistics as st, sys
sys.path.insert(0, 'validation')
from svgplot import Fig, Ax, COL

rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))
fl = float
pct = lambda v: f'{v * 100:.0f}%'
SYS = {'haploid': COL['blue'], 'diploid': COL['green'], 'haplodiploid': COL['amber']}
SYSM = {'haploid': 'o', 'diploid': 's', 'haplodiploid': '^'}
OUT = 'supplement_figures/'

# ---------- S2: power of SPOT against SLiM in each of the 20 runs, 6 scenarios with selection ----------
R = rd('validation/slim_power_runs.tsv')
scen = list(dict.fromkeys(r['scenario'] for r in R))
pal6 = [COL['blue'], COL['green'], COL['amber'], COL['purple'], COL['red'], COL['ink']]
name = lambda s: s.replace('_N', ', N = ').replace('_s', ', s = ').replace('_rec', ', recessive').replace('N = 1000', 'N = 1,000').replace('N = 2000', 'N = 2,000')
F = Fig(1000, 510)
a = Ax(F, 75, 50, 360, 360, (0.15, 0.8), (0.15, 0.8), title='Power in 20 independent runs per program', label=None)
a.frame([.2, .3, .4, .5, .6, .7, .8], [.2, .3, .4, .5, .6, .7, .8], 'power with SLiM', 'power with SPOT', xfmt=pct, yfmt=pct)
a.line([0.15, 0.8], [0.15, 0.8], COL['muted'], 1, '4 4')
for k, s in enumerate(scen):
    r = [x for x in R if x['scenario'] == s]
    a.points([fl(x['power_slim']) for x in r], [fl(x['power_spot']) for x in r], pal6[k], 'o', 3.2, True)
    mx, my = st.mean(fl(x['power_slim']) for x in r), st.mean(fl(x['power_spot']) for x in r)
    a.points([mx], [my], pal6[k], 'o', 5.5)
a.legend([(name(s), pal6[k], 'o') for k, s in enumerate(scen)] + [('hollow: one run; filled: mean of 20 runs', COL['muted'], 'oh')], a.x + a.w + 20, a.y + 18)
F.text(75, 490, 'Each run: threshold from 2,000 neutral simulations of the same program (tie rule of SPOT), power from 2,000 selection simulations; '
       'dashed: 1:1 (slim_compare_v2.js).', 11, color=COL['muted'])
print('FigS2', F.save(OUT + 'FigS2_slim_power_runs.svg'))
shutil.copy('validation/slim_regression.svg', OUT + 'FigS2b_slim_regression.svg')

# ---------- S3: SD of power between runs against the number of selection simulations ----------
M = rd('validation/mc_precision.tsv')
sc5 = list(dict.fromkeys(r['scenario'] for r in M))
F = Fig(1290, 450)
for j, (stat, title) in enumerate((('trend', 'Trend test'), ('two', 'First-versus-last test'))):
    ax = Ax(F, 75 + 470 * j, 50, 330, 300, (400, 6500), (0, 0.04), xlog=True, title=title, label='ab'[j])
    ax.frame([500, 1000, 2000, 5000], [0, .01, .02, .03, .04], 'selection simulations R', 'SD of power between 50 runs',
             xfmt=lambda v: f'{v:,.0f}', yfmt=lambda v: f'{v * 100:.0f} pt')
    for k, s in enumerate(sc5):
        r = sorted([x for x in M if x['scenario'] == s and x['statistic'] == stat], key=lambda x: fl(x['replicates']))
        xs = [fl(x['replicates']) for x in r]
        ax.line(xs, [fl(x['sd_power']) for x in r], pal6[k]); ax.points(xs, [fl(x['sd_power']) for x in r], pal6[k])
        ax.line(xs, [math.sqrt(fl(x['mean_power']) * (1 - fl(x['mean_power'])) / fl(x['replicates'])) for x in r], pal6[k], 1.2, '4 3')
    if j == 1:
        ax.legend([(s, pal6[k], 'o-') for k, s in enumerate(sc5)] + [('dashed: √(P(1 − P)/R), selection simulations only', COL['muted'], 'dash')],
                  ax.x + ax.w + 15, ax.y + 18)
F.text(75, 420, 'Example scenarios of the page; 50 runs with different seeds for each R, with 2R neutral simulations (monte_carlo_precision.js).', 11, color=COL['muted'])
print('FigS3', F.save(OUT + 'FigS3_precision.svg'))

# ---------- S4: rescaling, power with rescaling against full simulations; speed-up ----------
Rs = rd('validation/rescale_results.tsv')
F = Fig(1000, 470)
a = Ax(F, 75, 50, 330, 330, (0, 1), (0, 1), title='Power, trend test', label='a')
a.frame([0, .2, .4, .6, .8, 1], [0, .2, .4, .6, .8, 1], 'full simulation', 'with rescaling', xfmt=pct, yfmt=pct)
xs = [i / 100 for i in range(101)]
se = [2 * math.sqrt(2 * x * (1 - x) / 1000) for x in xs]
a.band(xs, [x - e for x, e in zip(xs, se)], [x + e for x, e in zip(xs, se)], COL['grey'], .25)
a.line([0, 1], [0, 1], COL['muted'], 1, '4 4')
for r in Rs:
    a.points([fl(r['power_exact'])], [fl(r['power_rescaled'])], SYS[r['system']], SYSM[r['system']], 3.4)
a.legend([(k, c, SYSM[k]) for k, c in SYS.items()], a.x + 10, a.y + 18)
Tm = rd('validation/rescale_times.tsv')   # run times measured alone (rescale_validation.js times), same order as Rs
assert [(r['G'], r['system'], r['N'], r['s'], r['p0']) for r in Tm] == [(r['G'], r['system'], r['N'], r['s'], r['p0']) for r in Rs]
sp = [fl(r['speedup']) for r in Tm]
b = Ax(F, 530, 50, 330, 330, (2, 40), (0.5, 100), xlog=True, ylog=True, title='Speed-up from rescaling', label='b')
b.frame([2, 5, 10, 20, 40], [1, 2, 5, 10, 20, 50, 100], 'rescaling factor λ', 'run time full / with rescaling')
b.line([2, 40], [2, 40], COL['muted'], 1, '4 4')
for r, v in zip(Rs, sp):
    b.points([fl(r['lambda'])], [v], SYS[r['system']], SYSM[r['system']], 3.4, True)
F.text(75, 430, f'{len(Rs)} test cases, 5,000–60,000 generations; grey band: ±2 Monte Carlo SE of a difference (1,000 selection simulations each); '
       f'b: dashed, speed-up = λ (rescale_validation.js).', 11, color=COL['muted'])
F.text(75, 446, f'Median speed-up {st.median(sp):.1f}, maximum {max(sp):.0f}; run times of one calculation each, measured with nothing else running (rescale_times.tsv).', 11, color=COL['muted'])
print('FigS4', F.save(OUT + 'FigS4_rescaling.svg'), 'median speed-up', round(st.median(sp), 2), 'max', round(max(sp), 1))

# ---------- S5: human ancient DNA (copy) ----------
shutil.copy('worked_examples/human_adna.svg', OUT + 'FigS5_human_adna.svg')

# ---------- S6: bacterial pathogen, power against s (rescaled) ----------
P = rd('worked_examples/pathogen.tsv')
F = Fig(840, 450)
f = Ax(F, 75, 50, 330, 300, (1.5e-5, 7e-4), (0, 1), xlog=True, title='Bacterial pathogen (17,000 generations, rescaled)', label=None)
f.frame([2e-5, 5e-5, 1e-4, 2e-4, 5e-4], [0, .2, .4, .6, .8, 1], 'selection coefficient s per generation', 'power', xfmt=lambda v: f'{v:g}', yfmt=pct)
f.hline(0.8, COL['muted'])
xs = [fl(x['s_per_generation']) for x in P]
f.line(xs, [fl(x['power_trend']) for x in P], COL['blue']); f.points(xs, [fl(x['power_trend']) for x in P], COL['blue'])
f.line(xs, [fl(x['power_first_vs_last']) for x in P], COL['green']); f.points(xs, [fl(x['power_first_vs_last']) for x in P], COL['green'], 's')
f.line(xs, [fl(x['power_naive']) for x in P], COL['red'], 1.2, '4 3')
f.legend([('trend test', COL['blue'], 'o-'), ('first-versus-last test', COL['green'], 's-'), ('sampling-only threshold (trend)', COL['red'], 'dash')], f.x + 10, f.y + 18)
F.text(75, 405, 'Haploid, 100 generations per year, N = 100,000, samples 1850/1900/1950/2000/2020 (10/10/15/20/40), p₀ = 0.05;', 11, color=COL['muted'])
F.text(75, 421, '500 selection simulations per point (worked_examples.js E).', 11, color=COL['muted'])
print('FigS6', F.save(OUT + 'FigS6_pathogen.svg'))
