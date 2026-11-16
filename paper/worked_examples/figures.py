######################################################################################
# Figures of the case studies (results from worked_examples.js and depth_case.js)    #
#                                                                                    #
# Three figures are drawn:                                                           #
#   bumblebees.svg : Figure 3 of the article (museum insect case study, 5 panels)    #
#   contrasts.svg  : long-lived mammal and bacterial pathogen (not in the article;   #
#                    the pathogen panel is drawn again as Figure S5 by               #
#                    paper/supplement_figures/make_figures.py)                       #
#   human_adna.svg : human ancient DNA (drawn again as Figure S6 by make_figures.py) #
#                                                                                    #
# Run from the repository root:                                                      #
#   python3 paper/worked_examples/figures.py                                         #
######################################################################################


#########################
# Load the dependencies #
#########################
import csv                                    # To read the tab-separated results files
import sys                                    # To add a folder to the Python search path
sys.path.insert(0, 'paper/validation')        # Folder that contains svgplot.py
from svgplot import Fig, Ax, COL, esc, STYLE  # Small plotting tools of this repository


###########################
# Small tools and colours #
###########################
rd = lambda f: list(csv.DictReader(open('paper/worked_examples/' + f), delimiter='\t'))   # Function: read a results file as a list of dictionaries (one per row)
fl = float                                    # Short name: text to number
pct = lambda v: f'{v * 100:.0f}%'             # Function: write a proportion as a percentage (tick labels)
pal = [COL['blue'], COL['green'], COL['amber'], COL['purple'], COL['red']]   # Colours used in turn for the curves of the two other figures
NECOL = ['#6B6B6B', '#000000']   # panel e: dark grey and black for the two Ne (blue and green are p0 colours)
P0COL = {'0.01': COL['blue'], '0.05': COL['green'], '0.1': COL['amber'], '0.3': COL['purple']}   # same colour for each p0 in every panel


####################################################################################
# Figure 3 of the article: museum insect case study (bumblebees.svg)               #
#                                                                                  #
# The legend of panels a to c is below these panels (colour = starting frequency;  #
# line style and marker = sample size). Panels d and e have their own legend on    #
# their right. There is no text under the figure: the settings are in the caption. #
####################################################################################
# ---------- Figure 1: bumblebees (Figure 3 of the paper) ----------
# Legend of a-c below these panels (colour = p0, line style + marker = sample size); panels d and e have their own legend on the right.
# No text under the figure: settings are in the caption of the paper.
TOP = 20                                   # panels start below the shared legend
DEFAULT = dict(STYLE)                         # Keep the default text sizes, to put them back after this figure
STYLE.update(tick=14, lab=17, title=15, title_centred=True)   # larger text for this figure, panel titles centred
LS = 14                                    # size of the legends and of the design names of panel d
F = Fig(1270, 1200)                           # Empty figure (width and height in pixels)
P3 = ['0.01', '0.1', '0.3']                   # The three starting frequencies shown
SC = {'1': ('o', None), '2': ('s', '7 3'), '5': ('^', '2 2.5')}   # sample-size factor: marker and line style
NE_P0 = [('0.01', '0.05'), ('0.1', '0.02'), ('0.3', '0.015')]       # (p0, s) of panels b-d

def key(ax, x, y, label, color, marker=None, dash=None, line=True, hollow=False, w=28):
    """one legend entry: optional line segment (with dash), optional marker, label"""
    if line:                                  # Short line segment
        da = f' stroke-dasharray="{dash}"' if dash else ''
        F.o.append(f'<line x1="{x}" x2="{x + w}" y1="{y - 4}" y2="{y - 4}" stroke="{color}" stroke-width="1.8"{da}/>')
    if marker:                                # Marker in the middle of the segment
        F.o.append(ax.marker(x + w / 2, y - 4, color, marker, 3.8, hollow))
    F.text(x + w + 7, y, esc(label), LS, color=COL['ink'])   # Text of the entry


####################################################
# Panel a: power against the selection coefficient #
####################################################
a = Ax(F, 90, 50 + TOP, 330, 290, (0.004, 0.25), (0, 1), xlog=True, title='Power against selection', label='a')   # Axes (x in log scale)

# legend of a-c, centred under the three panels (block about 860 px wide; panels span x = 90-1230)
LX = 230                                      # Left edge of the legend
F.text(LX, 446, 'Starting frequency (a–c):', LS, weight='bold', color=COL['ink'])   # First line of the legend: colours
for k, p in enumerate(P3):
    key(a, LX + 215 + k * 125, 446, f'p₀ = {p}', P0COL[p], 'o')
F.text(LX, 472, 'Individuals per date (b, c):', LS, weight='bold', color=COL['ink'])   # Second line: line styles and markers
for k, (sc, lab) in enumerate((('1', '20/20/20/40 (×1)'), ('2', '40/40/40/80 (×2)'), ('5', '100/100/100/200 (×5)'))):
    key(a, LX + 215 + k * 215, 472, lab, COL['ink'], SC[sc][0], SC[sc][1])

A = rd('bee_s_p0.tsv')                        # Results of part A of worked_examples.js
a.frame([0.005, 0.01, 0.02, 0.05, 0.1, 0.2], [0, .2, .4, .6, .8, 1], 'selection coefficient s', 'power (trend test)', yfmt=pct)   # Ticks and axis titles
a.hline(0.8, COL['muted'])                    # Dashed line at 80% power
for p0 in P3:                                 # One curve per starting frequency
    r = [x for x in A if x['p0'] == p0]
    xs = [fl(x['s']) for x in r]              # x values: selection coefficient
    ys = [fl(x['power_trend']) for x in r]    # y values: power of the trend test
    a.line(xs, ys, P0COL[p0])
    a.points(xs, ys, P0COL[p0])


#########################################################
# Panels b and c: power, and false-positive rate of the #
# sampling-only threshold, against Ne                   #
#########################################################
C = rd('bee_ne.tsv')                          # Results of part C of worked_examples.js
b = Ax(F, 535, 50 + TOP, 290, 290, (110, 2e5), (0, 1), xlog=True, title='Power against Ne', label='b')
b.frame([100, 1000, 10000, 100000], [0, .2, .4, .6, .8, 1], 'effective population size Ne', 'power (trend test)', yfmt=pct)
b.hline(0.8, COL['muted'])                    # Dashed line at 80% power
c = Ax(F, 940, 50 + TOP, 290, 290, (110, 2e5), (0, 0.4), xlog=True, title='False positives, sampling-only threshold', label='c')
c.frame([100, 1000, 10000, 100000], [0, .1, .2, .3, .4], 'effective population size Ne', 'false-positive rate (intended 5%)', yfmt=pct)
c.hline(0.05, COL['muted'])                   # Dashed line at the intended rate (5%)
for p0, s_ in NE_P0:                          # One colour per starting frequency
    for sc, (mk, dash) in SC.items():         # One line style and marker per sample size
        r = [x for x in C if x['p0'] == p0 and x['n_scale'] == sc]
        xs = [fl(x['Ne']) for x in r]         # x values: Ne
        for ax, col in ((b, 'power_trend'), (c, 'fpr_naive')):   # The same curves in panel b (power) and panel c (false positives)
            ys = [fl(x[col]) for x in r]
            ax.line(xs, ys, P0COL[p0], 1.6, dash)
            ax.points(xs, ys, P0COL[p0], mk, 3.4)


######################################################
# Panel d: gain in power from extra effort, compared #
# with the reference sampling                        #
######################################################
B = rd('bee_allocation.tsv')                  # Results of part B of worked_examples.js
designs = [x for x in list(dict.fromkeys(x['design'] for x in B))[1:] if '20x' not in x]  # historical 20x dropped: depth is in panel e
dlab = {'+60 historical: 40/40/40/40': '+60 historical (40/40/40/40)', '+60 modern: 20/20/20/100': '+60 modern (20/20/20/100)',
        'first and last only: 50/50 (same total)': 'first and last dates only (50/50)'}   # Names of the designs as written on the figure
d = Ax(F, 290, 525 + TOP, 560, 260, (-5, 25), (len(designs) - 0.5, -0.5), title='Gain in power from extra effort, compared with the reference sampling', label='d')   # Axes: one row per design
d.frame([-5, 0, 5, 10, 15, 20, 25], [], 'change in power (percentage points)', '', grid=False)
d.line([0, 0], [len(designs) - 0.5, -0.5], COL['muted'], 1.2, '5 4')   # Vertical dashed line at 0 (no change)
ref = {(x['Ne'], x['p0']): fl(x['power_trend']) for x in B if x['design'].startswith('reference')}   # Power of the reference sampling, for each Ne and p0
for i, lab in enumerate(designs):             # One row per design
    F.o.append(f'<line x1="{d.x}" x2="{d.x + d.w}" y1="{d.Y(i):.1f}" y2="{d.Y(i):.1f}" stroke="{COL["line"]}"/>')   # Thin horizontal guide line
    F.text(d.x - 8, d.Y(i) + 5, esc(dlab.get(lab, lab)), LS, 'end', COL['ink'])                                    # Name of the design, left of the panel
    for j, Ne in enumerate(['1500', '100000']):   # Strong drift (filled, slightly above the line), weak drift (hollow, slightly below)
        for p0 in P3:
            v = next(fl(x['power_trend']) for x in B if x['design'] == lab and x['Ne'] == Ne and x['p0'] == p0)   # Power of this design
            d.points([100 * (v - ref[(Ne, p0)])], [i + (-0.17 if j == 0 else 0.17)], P0COL[p0], 'o', 4.5, j == 1)   # Change in power, in percentage points

# Legend of panel d, on its right
dx, dy = d.x + d.w + 25, d.y + 18
F.text(dx, dy - 2, 'Panel d:', LS, weight='bold', color=COL['ink'])
for k, p in enumerate(P3):
    key(d, dx, dy + 23 + 23 * k, f'p₀ = {p}', P0COL[p], 'o', line=False)
key(d, dx, dy + 100, 'filled: Ne = 1,500 (strong drift)', COL['ink'], 'o', line=False)
key(d, dx, dy + 123, 'hollow: Ne = 100,000 (weak drift)', COL['ink'], 'o', line=False, hollow=True)


###########################################
# Panel e: power against sequencing depth #
###########################################
# (e) sequencing depth at all dates (depth_case.js), s = 0.02, p0 = 0.1, minimum depth 1
E_ = [x for x in rd('depth_case.tsv') if x['s'] == '0.02']   # Results of depth_case.js for s = 0.02
PERF = 60   # x position used for perfect genotypes
e = Ax(F, 290, 885 + TOP, 560, 230, (0.8, 80), (0, 1), xlog=True, title='Power against sequencing depth (p₀ = 0.1, s = 0.02)', label='e')
e.frame([1, 2, 5, 10, 30, PERF], [0, .2, .4, .6, .8, 1], 'mean depth at every sampling date', 'power (trend test)',
        xfmt=lambda v: 'perfect' if v == PERF else f'{v:g}×', yfmt=pct)   # Tick numbers written "1×", "2×", ..., and "perfect"
for k, Ne in enumerate(['1500', '100000']):   # One curve per level of drift
    col = NECOL[k]
    r = [x for x in E_ if x['Ne'] == Ne]
    get = lambda lab: next(fl(x['power']) for x in r if x['design'] == lab)   # Function: power of the design called "lab"
    xs = [1, 2, 5, 10, 30]                    # Depths
    ys = [get(f'{dd}x all dates') for dd in xs]
    e.line(xs, ys, col)
    e.points(xs, ys, col)
    e.points([PERF], [get('perfect')], col, 'd', 4.5)   # Perfect genotypes: diamond at the right of the panel
    for err, sh in (('0.02', 's'), ('0.05', '^')):      # Higher error rates per read, at 1x and 2x: hollow square and triangle
        e.points([1, 2], [get(f'1x all dates, error {err}'), get(f'2x all dates, error {err}')], col, sh, 4.2, True)

# Legend of panel e, on its right
ex, ey = e.x + e.w + 25, e.y + 18
F.text(ex, ey - 2, 'Panel e:', LS, weight='bold', color=COL['ink'])
key(e, ex, ey + 23, 'Ne = 1,500 (strong drift)', NECOL[0], 'o')
key(e, ex, ey + 46, 'Ne = 100,000 (weak drift)', NECOL[1], 'o')
key(e, ex, ey + 76, 'per-read error 0.005', COL['ink'], 'o', line=False)
key(e, ex, ey + 99, 'per-read error 0.02', COL['ink'], 's', line=False, hollow=True)
key(e, ex, ey + 122, 'per-read error 0.05', COL['ink'], '^', line=False, hollow=True)
key(e, ex, ey + 145, 'perfect genotypes', COL['ink'], 'd', line=False)

# Save the figure, then put the default text sizes back
print('bumblebees.svg', F.save('paper/worked_examples/bumblebees.svg'))
STYLE.clear()
STYLE.update(DEFAULT)   # the two figures below keep the default text sizes


############################################################
# Long-lived mammal and bacterial pathogen (contrasts.svg) #
############################################################
# ---------- Figure 2: contrasting examples ----------
F = Fig(1000, 440)                            # Empty figure

# Panel a: long-lived mammal, power against s for three demographic histories
M = rd('mammal.tsv')                          # Results of part D of worked_examples.js
e = Ax(F, 75, 50, 340, 290, (0, 0.52), (0, 1), title='Long-lived mammal (7 generations)', label='a')
e.frame([0, .1, .2, .3, .4, .5], [0, .2, .4, .6, .8, 1], 'selection coefficient s (per generation)', 'power', yfmt=pct)
e.hline(0.8, COL['muted'])                    # Dashed line at 80% power
hist = list(dict.fromkeys(x['history'] for x in M))   # The three histories, in the order of the file
order = [hist[1], hist[0], hist[2]]           # Drawn in this order: constant, bottleneck, severe bottleneck
for k, h in enumerate(order):                 # One colour per history
    r = [x for x in M if x['history'] == h]
    xs = [fl(x['s']) for x in r]
    e.line(xs, [fl(x['power_trend']) for x in r], pal[k])              # Solid: drift-aware threshold
    e.points(xs, [fl(x['power_trend']) for x in r], pal[k])
    e.line(xs, [fl(x['power_naive']) for x in r], pal[k], 1.2, '4 3')  # Dashed: sampling-only threshold
e.legend([(h.replace('bottleneck: ', 'bottleneck, ').replace('constant', 'constant'), pal[k], 'o-') for k, h in enumerate(order)]
         + [('dashed: sampling-only threshold', COL['muted'], 'dash')], e.x + 10, e.y + 18)

# Panel b: bacterial pathogen, power against s
P = rd('pathogen.tsv')                        # Results of part E of worked_examples.js
f = Ax(F, 570, 50, 340, 290, (1.5e-5, 7e-4), (0, 1), xlog=True, title='Bacterial pathogen (17,000 generations, rescaled)', label='b')
f.frame([2e-5, 5e-5, 1e-4, 2e-4, 5e-4], [0, .2, .4, .6, .8, 1], 'selection coefficient s (per generation)', 'power', yfmt=pct)
f.hline(0.8, COL['muted'])                    # Dashed line at 80% power
xs = [fl(x['s_per_generation']) for x in P]
f.line(xs, [fl(x['power_trend']) for x in P], pal[0])                  # Solid: drift-aware threshold
f.points(xs, [fl(x['power_trend']) for x in P], pal[0])
f.line(xs, [fl(x['power_naive']) for x in P], pal[0], 1.2, '4 3')      # Dashed: sampling-only threshold
f.legend([('drift-aware threshold', pal[0], 'o-'), ('sampling-only threshold', pal[0], 'dash')], f.x + 10, f.y + 18)

# Notes under the figure, then save
F.text(75, 392, 'a: diploid, generation time 20 years, samples 1880/1950/2020 (20/20/40), Ne 5,000 before 1950 and as shown after; p₀ = 0.1.', 11, color=COL['muted'])
F.text(75, 408, 'b: haploid, 100 generations per year, N = 100,000, samples 1850/1900/1950/2000/2020; p₀ = 0.05.', 11, color=COL['muted'])
F.text(75, 424, 'Example scenarios of the page; 1,000 (a) and 5,000 (b) replicates under selection per point.', 11, color=COL['muted'])
print('contrasts.svg', F.save('paper/worked_examples/contrasts.svg'))


######################################
# Human ancient DNA (human_adna.svg) #
######################################
# ---------- Figure 3: human ancient DNA ----------
F = Fig(1000, 440)                            # Empty figure

# Panel a: power against s, for several depths of the ancient samples
D = rd('human_depth.tsv')                     # Results of part F of worked_examples.js
g = Ax(F, 75, 50, 340, 290, (0.0008, 0.13), (0, 1), xlog=True, title='Human aDNA: depth of the ancient samples', label='a')
g.frame([0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1], [0, .2, .4, .6, .8, 1], 'selection coefficient s (per generation)', 'power (trend test)', yfmt=pct)
g.hline(0.8, COL['muted'])                    # Dashed line at 80% power
deps = ['1', '2', '5', '10', 'perfect']       # Depths of the ancient samples
for k, dp in enumerate(deps):                 # One colour per depth
    r = [x for x in D if x['ancient_depth'] == dp]
    xs = [fl(x['s']) for x in r]
    g.line(xs, [fl(x['power_trend']) for x in r], pal[k])
    g.points(xs, [fl(x['power_trend']) for x in r], pal[k])
g.legend([((dp + '×') if dp != 'perfect' else 'perfect genotypes', pal[k], 'o-') for k, dp in enumerate(deps)], g.x + 10, g.y + 18)

# Panel b: the same 100 reads per ancient date spread over more or fewer individuals
Bg = rd('human_budget.tsv')                   # Results of part G of worked_examples.js
h2 = Ax(F, 570, 50, 340, 290, (8, 125), (0, 1), xlog=True, title='Human aDNA: 100 reads per ancient sampling date', label='b')
h2.frame([10, 20, 25, 50, 100], [0, .2, .4, .6, .8, 1], 'individuals per ancient date (depth = 100 / individuals)', 'power (trend test)', yfmt=pct)
for k, s in enumerate(['0.005', '0.0075']):   # One colour per selection coefficient
    r = [x for x in Bg if x['s'] == s]
    xs = [fl(x['individuals_per_ancient_date']) for x in r]
    h2.line(xs, [fl(x['power_trend']) for x in r], pal[k])
    h2.points(xs, [fl(x['power_trend']) for x in r], pal[k])
h2.legend([(f's = {s}', pal[k], 'o-') for k, s in enumerate(['0.005', '0.0075'])], h2.x + 10, h2.y + 18)

# Notes under the figure, then save
F.text(75, 392, 'Diploid, generation time 28 years, samples at 7000, 5000, 3000 and 1000 BCE and 1000 and 2000 CE; Ne 10,000; p₀ = 0.1; present-day sample: 100 individuals at 20×.', 11, color=COL['muted'])
F.text(75, 408, 'a: 50 individuals per ancient date; b: the same 100 reads per ancient date spread over 100 individuals at 1× to 10 individuals at 10× (5,000 simulations with selection).', 11, color=COL['muted'])
F.text(75, 424, 'Drift-aware threshold; per-read error rate 0.01, minimum depth 1 read.', 11, color=COL['muted'])
print('human_adna.svg', F.save('paper/worked_examples/human_adna.svg'))
