#################################################################################
# Figure 2 of the article: false-positive rate of the two thresholds            #
#                                                                               #
# Three panels, all for the trend test:                                         #
#   a : false-positive rate against Ne, for the three inheritance systems       #
#   b : the reference design with one factor changed at a time                  #
#   c : many candidate loci (small alpha per locus)                             #
# The results come from false_positive_rate.js (fpr_A.tsv, fpr_B.tsv, fpr_C.tsv #
# and fpr_C_gaussian.tsv).                                                      #
#                                                                               #
# Run from the repository root:                                                 #
#   python3 paper/validation/fpr_figure.py                                      #
# It writes paper/validation/false_positive_rate.svg.                           #
#################################################################################


#########################
# Load the dependencies #
#########################
import csv                                    # To read the tab-separated results files
import sys                                    # To add a folder to the Python search path
sys.path.insert(0, 'paper/validation')        # Folder that contains svgplot.py
from svgplot import Fig, Ax, COL, esc, STYLE  # Small plotting tools of this repository
STYLE.update(tick=14, lab=17)   # larger tick numbers and axis titles than the default (11 and 12)


#########################################
# Read the results and start the figure #
#########################################
rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))   # Function: read a tab-separated file as a list of dictionaries (one per row)
A_ = rd('paper/validation/fpr_A.tsv')         # Part A: Ne x inheritance system
B_ = rd('paper/validation/fpr_B.tsv')         # Part B: one design factor at a time
C_ = rd('paper/validation/fpr_C.tsv')         # Part C: many candidate loci
G_ = rd('paper/validation/fpr_C_gaussian.tsv')   # same analysis with a normal fit to 2,000 neutral simulations, for comparison
fl = lambda r, k: float(r[k])                 # Function: value of column k in row r, as a number
F = Fig(1285, 440)                            # Empty figure (width and height in pixels)
SH = {'haploid': 'o', 'diploid': 's', 'haplodiploid': '^'}   # Marker shape of each inheritance system (circle, square, triangle)


####################################################################################
# Panel a: false-positive rate against Ne                                          #
# Red: sampling-only threshold; blue: drift-aware threshold (mean of 20 runs, with #
# a bar of +/- 1 standard deviation)                                               #
####################################################################################
# (a) false-positive rate against Ne, trend test
a = Ax(F, 85, 50, 310, 320, (70, 1.4e6), (0, 0.2), xlog=True, title='False-positive rate against Ne', label='a')   # Axes: position, size, limits of x (log scale) and y
a.frame([100, 1000, 10000, 100000, 1000000], [0, .05, .1, .15, .2], 'effective population size Ne', 'false-positive rate (α = 0.05)')   # Ticks and axis titles
a.hline(0.05, COL['muted'])                   # Dashed line at the intended rate (5%)
for sysm, sh in SH.items():                   # One pair of curves per inheritance system
    rows = [r for r in A_ if r['system'] == sysm and r['statistic'] == 'trend']   # Results of this system, trend test
    ne = [fl(r, 'Ne') for r in rows]          # x values: Ne
    a.line(ne, [fl(r, 'fpr_naive') for r in rows], COL['red'], 1.4)               # Sampling-only threshold: line
    a.points(ne, [fl(r, 'fpr_naive') for r in rows], COL['red'], sh, 4)           # Sampling-only threshold: markers
    m = [fl(r, 'fpr_drift_aware_mean') for r in rows]                             # Drift-aware threshold: mean over the 20 runs
    s = [fl(r, 'fpr_drift_aware_sd') for r in rows]                               # Drift-aware threshold: standard deviation between runs
    a.errbar(ne, [x - y for x, y in zip(m, s)], [x + y for x, y in zip(m, s)], COL['blue'])   # Bars of +/- 1 standard deviation
    a.line(ne, m, COL['blue'], 1.4)
    a.points(ne, m, COL['blue'], sh, 4)
a.legend([('sampling-only threshold', COL['red'], 'o-'), ('drift-aware threshold (±1 SD)', COL['blue'], 'o-'),
          ('haploid', COL['muted'], 'o'), ('diploid', COL['muted'], 's'), ('haplodiploid', COL['muted'], '^')], a.X(4000), a.Y(0.175))   # Legend, placed inside the panel


##################################################################################
# Panel b: the reference design with one factor changed at a time                #
# One line per design; blue: drift-aware threshold; red: sampling-only threshold #
##################################################################################
# (b) one design factor at a time (haplodiploid, Ne = 1,500), trend test
rows = [r for r in B_ if r['statistic'] == 'trend']   # Results of the trend test
lab = lambda r: r['label'][2:].replace('base: 4 points, 20/20/20/40, 10x/15x, p0 = 0.1', 'reference design').replace('x', '×').replace('p0', 'p₀')   # Function: name of the design as written on the figure
b = Ax(F, 630, 50, 230, 320, (0, 0.25), (len(rows) - 0.5, -0.5), title='Design variations (Ne = 1,500)', label='b')   # Axes: y goes from the first design (top) to the last (bottom)
b.frame([0, .05, .1, .15, .2, .25], [], 'false-positive rate (α = 0.05)', '', grid=False)   # Ticks of x only
b.line([0.05, 0.05], [len(rows) - 0.5, -0.5], COL['muted'], 1.2, '5 4')                     # Vertical dashed line at the intended rate (5%)
for i, r in enumerate(rows):                  # One line per design
    y = i
    F.o.append(f'<line x1="{b.x}" x2="{b.x + b.w}" y1="{b.Y(y):.1f}" y2="{b.Y(y):.1f}" stroke="{COL["line"]}"/>')   # Thin horizontal guide line
    F.text(b.x - 8, b.Y(y) + 4, esc(lab(r)), 13, 'end', COL['ink'])                                                  # Name of the design, left of the panel
    m = fl(r, 'fpr_drift_aware_mean')         # Drift-aware threshold: mean over the 20 runs
    s = fl(r, 'fpr_drift_aware_sd')           # Drift-aware threshold: standard deviation between runs
    F.o.append(f'<line x1="{b.X(m - s):.1f}" x2="{b.X(m + s):.1f}" y1="{b.Y(y):.1f}" y2="{b.Y(y):.1f}" stroke="{COL["blue"]}" stroke-width="1.4"/>')   # Bar of +/- 1 standard deviation
    b.points([m], [y], COL['blue'], 'o', 4)                       # Drift-aware threshold
    b.points([fl(r, 'fpr_naive')], [y], COL['red'], 'o', 4)       # Sampling-only threshold


#########################################################################################
# Panel c: many candidate loci                                                          #
# y: realised false-positive rate divided by the intended rate (1 = as intended)        #
# Grey dashed: threshold from a normal distribution fitted to 2,000 neutral simulations #
#########################################################################################
# (c) small per-locus alpha: realised / nominal false-positive rate, trend test
c = Ax(F, 975, 50, 280, 320, (3e-5, 0.08), (0.1, 1000), xlog=True, ylog=True, title='Many loci tested (Bonferroni)', label='c')   # Axes: both in log scale
import math
pw = lambda v: f'5×10<tspan dy="-6" font-size="10">−{-int(math.floor(math.log10(v)))}</tspan>'   # 5×10⁻⁵ ... 5×10⁻²
c.frame([5e-5, 5e-4, 5e-3, 5e-2], [0.1, 1, 10, 100, 1000], 'per-locus α′ = 0.05 / L', 'realised / intended false-positive rate', xfmt=pw)   # Ticks and axis titles
F.o.append(f'<rect x="{c.X(3e-5):.1f}" y="{c.y}" width="{c.X(0.0025) - c.X(3e-5):.1f}" height="{c.h}" fill="{COL["grey"]}" fill-opacity=".12"/>')   # Grey area: SPOT runs more neutral simulations there
F.text(c.X(3.5e-5), c.Y(0.25), 'more neutral simulations', 10, color=COL['muted'])   # Note inside the grey area (two lines)
F.text(c.X(3.5e-5), c.Y(0.25) + 12, '(20,000 and 200,000)', 10, color=COL['muted'])
c.hline(1, COL['muted'])                      # Dashed line at 1 (realised rate = intended rate)
DS = [('haplodiploid Ne=300', 'o', 'Ne = 300'), ('haplodiploid Ne=1500', 's', 'Ne = 1,500'), ('haplodiploid Ne=100000', '^', 'Ne = 100,000'),
      ('haplodiploid Ne=1500, p0=0.01', 'd', 'Ne = 1,500, p₀ 0.01')]   # The four designs: (name in the results file, marker shape, name in the legend)
for key, sh, _ in DS:                         # One set of curves per design
    rows = [r for r in C_ if r['label'].startswith('C ' + key + ' m=') and r['statistic'] == 'trend']   # Results of part C for this design
    al = [fl(r, 'alpha_per_locus') for r in rows]                                                      # x values: alpha per locus
    old = [r for r in G_ if r['label'].startswith('G ' + key + ' m=') and r['statistic'] == 'trend']   # Results with the normal fit to 2,000 neutral simulations
    ys = [max(0.1, fl(r, 'fpr_drift_aware_mean') / x) for r, x in zip(old, al)]   # Realised / intended rate (values below 0.1 are drawn at 0.1, the bottom of the axis)
    c.line(al, ys, COL['grey'], 1.1, '4 3')
    c.points(al, ys, COL['grey'], sh, 3.4, True)
    for col, k in ((COL['red'], 'fpr_naive'), (COL['blue'], 'fpr_drift_aware_mean')):   # Red: sampling-only threshold; blue: drift-aware threshold
        ys = [max(0.1, fl(r, k) / x) for r, x in zip(rows, al)]
        c.line(al, ys, col, 1.2)
        c.points(al, ys, col, sh, 3.8)
c.legend([(lab_, COL['muted'], sh) for _, sh, lab_ in DS] + [('normal fit (2,000 neutral)', COL['grey'], 'dash')], c.X(0.0011), c.Y(700))   # Legend


###################
# Save the figure #
###################
# no text under the figure: the settings are in the caption of the paper
print('figure bytes', F.save('paper/validation/false_positive_rate.svg'))   # Save and print the size of the file
