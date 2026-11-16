################################################################################################
# Regression of the SLiM results on the SPOT results (Supporting Information S4, Figure S2)    #
#                                                                                              #
# For three quantities (mean allele frequency, power, false-positive rate), the value given    #
# by SLiM is regressed on the value given by SPOT (ordinary least squares).                    #
# If the two programs agree, the slope is 1 and the intercept is 0.                            #
#                                                                                              #
# Inputs:                                                                                      #
#   paper/validation/slim_results.tsv    (made by slim_compare.js)                             #
#   paper/validation/slim_power_v2.tsv   (made by slim_compare_v2.js)                          #
#   paper/validation/slim_fpr_v2.tsv     (made by slim_compare_v2.js)                          #
# Outputs:                                                                                     #
#   paper/validation/slim_regression.tsv (slope, intercept, R2, with 95% confidence intervals) #
#   paper/validation/slim_regression.svg (the three scatter plots)                             #
#                                                                                              #
# Run from the repository root:                                                                #
#   python3 paper/validation/slim_regression.py                                                #
################################################################################################


#########################
# Load the dependencies #
#########################
import csv                             # To read the tab-separated results files
import sys                             # To add a folder to the Python search path
import math                            # For the square root
sys.path.insert(0, 'paper/validation')   # Folder that contains svgplot.py
from svgplot import Fig, Ax, COL       # Small plotting tools of this repository (figure, axes, colours)


################################
# Small tools and fixed values #
################################
def rd(f):
    """Read a tab-separated file and return its rows as a list of dictionaries (one per row)."""
    return list(csv.DictReader(open(f), delimiter='\t'))

def sysof(name):
    """Return the inheritance system from a scenario name (for example 'diploid_N1000_s0.02' gives 'diploid')."""
    return name.split('_')[0]

SH = {'haploid': 'o', 'diploid': 's', 'haplodiploid': '^'}   # Marker shape of each inheritance system (circle, square, triangle)
T975 = {8: 2.306, 31: 2.040, 4: 2.776}                       # 97.5% quantile of Student's t, for n - 2 degrees of freedom (95% confidence intervals)


#################################################################
# Function: ordinary least squares regression of y on x         #
# Returns the slope and the intercept with their 95% confidence #
# intervals, R2 and the largest absolute residual               #
#################################################################
def ols(x, y):
    n = len(x)                                               # Number of points
    mx = sum(x) / n                                          # Mean of x
    my = sum(y) / n                                          # Mean of y
    sxx = sum((a - mx) ** 2 for a in x)                      # Sum of squares of x
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))     # Sum of cross-products of x and y
    b = sxy / sxx                                            # Slope
    a = my - b * mx                                          # Intercept
    res = [yy - a - b * xx for xx, yy in zip(x, y)]          # Residuals (observed y minus fitted y)
    sse = sum(r * r for r in res)                            # Sum of squared residuals
    syy = sum((v - my) ** 2 for v in y)                      # Sum of squares of y
    s2 = sse / (n - 2)                                       # Residual variance
    t = T975[n - 2]                                          # Student's t quantile for this number of points
    seb = math.sqrt(s2 / sxx)                                # Standard error of the slope
    sea = math.sqrt(s2 * (1 / n + mx * mx / sxx))            # Standard error of the intercept
    return dict(n=n,
                slope=b,
                slope_lo=b - t * seb,                        # Lower limit of the 95% confidence interval of the slope
                slope_hi=b + t * seb,                        # Upper limit
                intercept=a,
                int_lo=a - t * sea,                          # Lower limit of the 95% confidence interval of the intercept
                int_hi=a + t * sea,                          # Upper limit
                r2=1 - sse / syy,                            # Proportion of the variance of y explained by x
                max_abs_residual=max(abs(r) for r in res))   # Largest distance between a point and the fitted line


##############################################################
# Read the results of the two programs                       #
# Each point is (SPOT value, SLiM value, inheritance system, #
# group used to choose the colour)                           #
##############################################################
# Panel a: mean allele frequency at generations 20, 40 and 90
F_ = rd('paper/validation/slim_results.tsv')
fa = [(float(r['mean_spot']), float(r['mean_slim']), sysof(r['scenario']), r['run']) for r in F_]

# Panel b: mean power over 20 runs, for the scenarios with selection
P_ = rd('paper/validation/slim_power_v2.tsv')
pb = [(float(r['mean_power_spot']), float(r['mean_power_slim']), r['system'], 'selected') for r in P_]

# Panel c: false-positive rate of the two thresholds, for the neutral scenarios
R_ = rd('paper/validation/slim_fpr_v2.tsv')
fc_spot = [(float(r['fpr_spot_threshold_on_spot']), float(r['fpr_spot_threshold_on_slim']), r['system'], 'spot') for r in R_]    # Drift-aware threshold
fc_naive = [(float(r['fpr_sampling_only_on_spot']), float(r['fpr_sampling_only_on_slim']), r['system'], 'naive') for r in R_]   # Sampling-only threshold
fc = fc_spot + fc_naive


####################################################
# One regression and one scatter plot per quantity #
####################################################
out = [['panel', 'quantity', 'n', 'slope', 'slope_95lo', 'slope_95hi', 'intercept', 'intercept_95lo', 'intercept_95hi', 'R2', 'max_abs_residual']]   # Header of the results table
F = Fig(1060, 410)                                           # Empty figure (width and height in pixels)

# For each panel: data, axis limits, axis ticks, title, panel letter, description written in the results table
panels = [
    (fa, (0, 0.5), [0, .1, .2, .3, .4, .5], 'Mean allele frequency', 'a', 'mean frequency at generations 20, 40, 90 (11 scenarios)'),
    (pb, (0, 1), [0, .2, .4, .6, .8, 1], 'Power', 'b', 'mean power over 20 runs (6 scenarios with selection)'),
    (fc, (0, 0.25), [0, .05, .1, .15, .2, .25], 'False-positive rate', 'c', 'false-positive rate, both thresholds (5 neutral scenarios)')]

for k, (data, lim, ticks, title, lab, q) in enumerate(panels):
    x = [d[0] for d in data]                                 # SPOT values
    y = [d[1] for d in data]                                 # SLiM values
    r = ols(x, y)                                            # Regression of SLiM on SPOT

    # Add one line to the results table
    out.append([lab, q, r['n']] + [f"{r[c]:.4f}" for c in ('slope', 'slope_lo', 'slope_hi', 'intercept', 'int_lo', 'int_hi', 'r2', 'max_abs_residual')])

    # Draw the panel: axes, 1:1 line, points, slope and R2
    ax = Ax(F, 75 + k * 340, 50, 250, 250, lim, lim, title=title, label=lab)   # Axes of panel k (the panels are side by side)
    if lab != 'a':
        fmt = lambda v: f'{v:.0%}'                           # Panels b and c: tick labels as percentages
    else:
        fmt = lambda v: f'{v:g}'                             # Panel a: tick labels as plain numbers
    ax.frame(ticks, ticks, 'SPOT', 'SLiM', xfmt=fmt, yfmt=fmt)                 # Frame, grid and axis titles
    ax.line([lim[0], lim[1]], [lim[0], lim[1]], COL['muted'], 1.2, '5 4')      # Dashed 1:1 line (perfect agreement)
    for xx, yy, s_, g in data:                               # One marker per point
        col = {'selected': COL['blue'], 'neutral': COL['muted'], 'spot': COL['blue'], 'naive': COL['red']}[g]   # Colour of the group
        ax.points([xx], [yy], col, SH[s_], 4, g == 'neutral')                  # Shape = inheritance system; hollow if neutral
    F.text(ax.x + 8, ax.y + 18, f"slope {r['slope']:.3f}, R² {r['r2']:.3f}", 11, 'start', COL['ink'])   # Slope and R2 written in the panel


##################################
# Notes written under the figure #
##################################
F.text(75, 372, 'Dashed line: 1:1. Shapes: circle haploid, square diploid, triangle haplodiploid. (a) filled, with selection; hollow, neutral. '
       '(b) each program\'s own threshold.', 11, color=COL['muted'])
F.text(75, 388, '(c) blue, threshold from 20,000 SPOT neutral simulations; red, sampling-only threshold; each applied to 10,000 neutral simulations of each program.', 11, color=COL['muted'])


#################################
# Save the figure and the table #
#################################
print('slim_regression.svg', F.save('paper/validation/slim_regression.svg'))   # Save the figure and print its size
with open('paper/validation/slim_regression.tsv', 'w') as f:                   # Save the results table (tab-separated)
    f.write('\n'.join('\t'.join(map(str, r)) for r in out) + '\n')
for r in out:                                                                  # Also print the table on the screen
    print('\t'.join(map(str, r)))
