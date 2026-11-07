"""Figures of the worked examples (results from worked_examples.js).
Run from the repository root: python3 paper/worked_examples/figures.py  ->  paper/worked_examples/bumblebees.svg, contrasts.svg"""
import csv, sys
sys.path.insert(0, 'paper/validation')
from svgplot import Fig, Ax, COL, esc

rd = lambda f: list(csv.DictReader(open('paper/worked_examples/' + f), delimiter='\t'))
fl = float
pct = lambda v: f'{v * 100:.0f}%'
pal = [COL['blue'], COL['green'], COL['amber'], COL['purple'], COL['red']]
NECOL = ['#6B6B6B', '#000000']   # panel e: dark grey and black for the two Ne (blue and green are p0 colours)
P0COL = {'0.01': COL['blue'], '0.05': COL['green'], '0.1': COL['amber'], '0.3': COL['purple']}   # same colour for each p0 in every panel

# ---------- Figure 1: bumblebees (Figure 3 of the paper) ----------
# Legend of a-c below these panels (colour = p0, line style + marker = sample size); panels d and e have their own legend on the right.
# No text under the figure: settings are in the caption of the paper.
TOP = 20                                   # panels start below the shared legend
F = Fig(1250, 1140)
P3 = ['0.01', '0.1', '0.3']
SC = {'1': ('o', None), '2': ('s', '7 3'), '5': ('^', '2 2.5')}   # sample-size factor: marker and line style
NE_P0 = [('0.01', '0.05'), ('0.1', '0.02'), ('0.3', '0.015')]       # (p0, s) of panels b-d

def key(ax, x, y, label, color, marker=None, dash=None, line=True, hollow=False, w=28):
    """one legend entry: optional line segment (with dash), optional marker, label"""
    if line:
        da = f' stroke-dasharray="{dash}"' if dash else ''
        F.o.append(f'<line x1="{x}" x2="{x + w}" y1="{y - 4}" y2="{y - 4}" stroke="{color}" stroke-width="1.8"{da}/>')
    if marker:
        F.o.append(ax.marker(x + w / 2, y - 4, color, marker, 3.8, hollow))
    F.text(x + w + 7, y, esc(label), 12, color=COL['ink'])

a = Ax(F, 75, 50 + TOP, 330, 290, (0.004, 0.25), (0, 1), xlog=True, title='Power against selection', label='a')
# legend of a-c, centred under the three panels (block about 745 px wide; panels span x = 75-1210)
LX = 195
F.text(75 + LX, 432, 'Starting frequency (a–c):', 12, weight='bold', color=COL['ink'])
for k, p in enumerate(P3):
    key(a, 290 + LX + k * 105, 432, f'p₀ = {p}', P0COL[p], 'o')
F.text(75 + LX, 456, 'Individuals per date (b, c):', 12, weight='bold', color=COL['ink'])
for k, (sc, lab) in enumerate((('1', '20/20/20/40 (×1)'), ('2', '40/40/40/80 (×2)'), ('5', '100/100/100/200 (×5)'))):
    key(a, 290 + LX + k * 190, 456, lab, COL['ink'], SC[sc][0], SC[sc][1])

A = rd('bee_s_p0.tsv')
a.frame([0.005, 0.01, 0.02, 0.05, 0.1, 0.2], [0, .2, .4, .6, .8, 1], 'selection coefficient s', 'power (trend test)', yfmt=pct)
a.hline(0.8, COL['muted'])
for p0 in P3:
    r = [x for x in A if x['p0'] == p0]
    xs, ys = [fl(x['s']) for x in r], [fl(x['power_trend']) for x in r]
    a.line(xs, ys, P0COL[p0]); a.points(xs, ys, P0COL[p0])

C = rd('bee_ne.tsv')
b = Ax(F, 520, 50 + TOP, 290, 290, (110, 2e5), (0, 1), xlog=True, title='Power against Ne', label='b')
b.frame([100, 1000, 10000, 100000], [0, .2, .4, .6, .8, 1], 'effective population size Ne', 'power (trend test)', yfmt=pct)
b.hline(0.8, COL['muted'])
c = Ax(F, 920, 50 + TOP, 290, 290, (110, 2e5), (0, 0.4), xlog=True, title='False positives, sampling-only threshold', label='c')
c.frame([100, 1000, 10000, 100000], [0, .1, .2, .3, .4], 'effective population size Ne', 'false-positive rate (intended 5%)', yfmt=pct)
c.hline(0.05, COL['muted'])
for p0, s_ in NE_P0:
    for sc, (mk, dash) in SC.items():
        r = [x for x in C if x['p0'] == p0 and x['n_scale'] == sc]
        xs = [fl(x['Ne']) for x in r]
        for ax, col in ((b, 'power_trend'), (c, 'fpr_naive')):
            ys = [fl(x[col]) for x in r]
            ax.line(xs, ys, P0COL[p0], 1.6, dash); ax.points(xs, ys, P0COL[p0], mk, 3.4)

B = rd('bee_allocation.tsv')
designs = [x for x in list(dict.fromkeys(x['design'] for x in B))[1:] if '20x' not in x]  # historical 20x dropped: depth is in panel e
dlab = {'+60 historical: 40/40/40/40': '+60 historical (40/40/40/40)', '+60 modern: 20/20/20/100': '+60 modern (20/20/20/100)',
        'first and last only: 50/50 (same total)': 'first and last dates only (50/50)'}
d = Ax(F, 260, 495 + TOP, 560, 260, (-5, 25), (len(designs) - 0.5, -0.5), title='Gain in power from extra effort, compared with the reference sampling', label='d')
d.frame([-5, 0, 5, 10, 15, 20, 25], [], 'change in power (percentage points)', '', grid=False)
d.line([0, 0], [len(designs) - 0.5, -0.5], COL['muted'], 1.2, '5 4')
ref = {(x['Ne'], x['p0']): fl(x['power_trend']) for x in B if x['design'].startswith('reference')}
for i, lab in enumerate(designs):
    F.o.append(f'<line x1="{d.x}" x2="{d.x + d.w}" y1="{d.Y(i):.1f}" y2="{d.Y(i):.1f}" stroke="{COL["line"]}"/>')
    F.text(d.x - 8, d.Y(i) + 4, esc(dlab.get(lab, lab)), 12, 'end', COL['ink'])
    for j, Ne in enumerate(['1500', '100000']):
        for p0 in P3:
            v = next(fl(x['power_trend']) for x in B if x['design'] == lab and x['Ne'] == Ne and x['p0'] == p0)
            d.points([100 * (v - ref[(Ne, p0)])], [i + (-0.17 if j == 0 else 0.17)], P0COL[p0], 'o', 4.5, j == 1)

dx, dy = d.x + d.w + 25, d.y + 18
F.text(dx, dy - 2, 'Panel d:', 12, weight='bold', color=COL['ink'])
for k, p in enumerate(P3):
    key(d, dx, dy + 20 + 20 * k, f'p₀ = {p}', P0COL[p], 'o', line=False)
key(d, dx, dy + 86, 'filled: Ne = 1,500 (strong drift)', COL['ink'], 'o', line=False)
key(d, dx, dy + 106, 'hollow: Ne = 100,000 (weak drift)', COL['ink'], 'o', line=False, hollow=True)

# (e) sequencing depth at all dates (depth_case.js), s = 0.02, p0 = 0.1, minimum depth 1
E_ = [x for x in rd('depth_case.tsv') if x['s'] == '0.02']
PERF = 60   # x position used for perfect genotypes
e = Ax(F, 260, 840 + TOP, 560, 230, (0.8, 80), (0, 1), xlog=True, title='Power against sequencing depth (p₀ = 0.1, s = 0.02)', label='e')
e.frame([1, 2, 5, 10, 30, PERF], [0, .2, .4, .6, .8, 1], 'mean depth at every sampling date', 'power (trend test)',
        xfmt=lambda v: 'perfect' if v == PERF else f'{v:g}×', yfmt=pct)
for k, Ne in enumerate(['1500', '100000']):
    col = NECOL[k]
    r = [x for x in E_ if x['Ne'] == Ne]
    get = lambda lab: next(fl(x['power']) for x in r if x['design'] == lab)
    xs = [1, 2, 5, 10, 30]
    ys = [get(f'{dd}x all dates') for dd in xs]
    e.line(xs, ys, col); e.points(xs, ys, col); e.points([PERF], [get('perfect')], col, 'd', 4.5)
    for err, sh in (('0.02', 's'), ('0.05', '^')):
        e.points([1, 2], [get(f'1x all dates, error {err}'), get(f'2x all dates, error {err}')], col, sh, 4.2, True)
ex, ey = e.x + e.w + 25, e.y + 18
F.text(ex, ey - 2, 'Panel e:', 12, weight='bold', color=COL['ink'])
key(e, ex, ey + 20, 'Ne = 1,500 (strong drift)', NECOL[0], 'o')
key(e, ex, ey + 40, 'Ne = 100,000 (weak drift)', NECOL[1], 'o')
key(e, ex, ey + 66, 'per-read error 0.005', COL['ink'], 'o', line=False)
key(e, ex, ey + 86, 'per-read error 0.02', COL['ink'], 's', line=False, hollow=True)
key(e, ex, ey + 106, 'per-read error 0.05', COL['ink'], '^', line=False, hollow=True)
key(e, ex, ey + 126, 'perfect genotypes', COL['ink'], 'd', line=False)
print('bumblebees.svg', F.save('paper/worked_examples/bumblebees.svg'))

# ---------- Figure 2: contrasting examples ----------
F = Fig(1000, 440)
M = rd('mammal.tsv')
e = Ax(F, 75, 50, 340, 290, (0, 0.52), (0, 1), title='Long-lived mammal (7 generations)', label='a')
e.frame([0, .1, .2, .3, .4, .5], [0, .2, .4, .6, .8, 1], 'selection coefficient s (per generation)', 'power', yfmt=pct)
e.hline(0.8, COL['muted'])
hist = list(dict.fromkeys(x['history'] for x in M))
order = [hist[1], hist[0], hist[2]]
for k, h in enumerate(order):
    r = [x for x in M if x['history'] == h]
    xs = [fl(x['s']) for x in r]
    e.line(xs, [fl(x['power_trend']) for x in r], pal[k]); e.points(xs, [fl(x['power_trend']) for x in r], pal[k])
    e.line(xs, [fl(x['power_naive']) for x in r], pal[k], 1.2, '4 3')
e.legend([(h.replace('bottleneck: ', 'bottleneck, ').replace('constant', 'constant'), pal[k], 'o-') for k, h in enumerate(order)]
         + [('dashed: sampling-only threshold', COL['muted'], 'dash')], e.x + 10, e.y + 18)
P = rd('pathogen.tsv')
f = Ax(F, 570, 50, 340, 290, (1.5e-5, 7e-4), (0, 1), xlog=True, title='Bacterial pathogen (17,000 generations, rescaled)', label='b')
f.frame([2e-5, 5e-5, 1e-4, 2e-4, 5e-4], [0, .2, .4, .6, .8, 1], 'selection coefficient s (per generation)', 'power', yfmt=pct)
f.hline(0.8, COL['muted'])
xs = [fl(x['s_per_generation']) for x in P]
f.line(xs, [fl(x['power_trend']) for x in P], pal[0]); f.points(xs, [fl(x['power_trend']) for x in P], pal[0])
f.line(xs, [fl(x['power_naive']) for x in P], pal[0], 1.2, '4 3')
f.legend([('drift-aware threshold', pal[0], 'o-'), ('sampling-only threshold', pal[0], 'dash')], f.x + 10, f.y + 18)
F.text(75, 392, 'a: diploid, generation time 20 years, samples 1880/1950/2020 (20/20/40), Ne 5,000 before 1950 and as shown after; p₀ = 0.1.', 11, color=COL['muted'])
F.text(75, 408, 'b: haploid, 100 generations per year, N = 100,000, samples 1850/1900/1950/2000/2020; p₀ = 0.05.', 11, color=COL['muted'])
F.text(75, 424, 'Example scenarios of the page; 1,000 (a) and 500 (b) replicates under selection per point.', 11, color=COL['muted'])
print('contrasts.svg', F.save('paper/worked_examples/contrasts.svg'))

# ---------- Figure 3: human ancient DNA ----------
F = Fig(1000, 440)
D = rd('human_depth.tsv')
g = Ax(F, 75, 50, 340, 290, (0.0008, 0.13), (0, 1), xlog=True, title='Human aDNA: depth of the ancient samples', label='a')
g.frame([0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1], [0, .2, .4, .6, .8, 1], 'selection coefficient s (per generation)', 'power (trend test)', yfmt=pct)
g.hline(0.8, COL['muted'])
deps = ['1', '2', '5', '10', 'perfect']
for k, dp in enumerate(deps):
    r = [x for x in D if x['ancient_depth'] == dp]
    xs = [fl(x['s']) for x in r]
    g.line(xs, [fl(x['power_trend']) for x in r], pal[k]); g.points(xs, [fl(x['power_trend']) for x in r], pal[k])
g.legend([((dp + '×') if dp != 'perfect' else 'perfect genotypes', pal[k], 'o-') for k, dp in enumerate(deps)], g.x + 10, g.y + 18)
Bg = rd('human_budget.tsv')
h2 = Ax(F, 570, 50, 340, 290, (8, 125), (0, 1), xlog=True, title='Human aDNA: 100 reads per ancient sampling date', label='b')
h2.frame([10, 20, 25, 50, 100], [0, .2, .4, .6, .8, 1], 'individuals per ancient date (depth = 100 / individuals)', 'power (trend test)', yfmt=pct)
for k, s in enumerate(['0.005', '0.0075']):
    r = [x for x in Bg if x['s'] == s]
    xs = [fl(x['individuals_per_ancient_date']) for x in r]
    h2.line(xs, [fl(x['power_trend']) for x in r], pal[k]); h2.points(xs, [fl(x['power_trend']) for x in r], pal[k])
h2.legend([(f's = {s}', pal[k], 'o-') for k, s in enumerate(['0.005', '0.0075'])], h2.x + 10, h2.y + 18)
F.text(75, 392, 'Diploid, generation time 28 years, samples at 7000, 5000, 3000 and 1000 BCE and 1000 and 2000 CE; Ne 10,000; p₀ = 0.1; present-day sample: 100 individuals at 20×.', 11, color=COL['muted'])
F.text(75, 408, 'a: 50 individuals per ancient date; b: the same 100 reads per ancient date spread over 100 individuals at 1× to 10 individuals at 10× (5,000 simulations with selection).', 11, color=COL['muted'])
F.text(75, 424, 'Drift-aware threshold; per-read error rate 0.01, minimum depth 1 read.', 11, color=COL['muted'])
print('human_adna.svg', F.save('paper/worked_examples/human_adna.svg'))
