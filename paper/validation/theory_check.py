##########################################################################################
# Agreement of SPOT with theory (Supporting Information S3, Table S4)                    #
#                                                                                        #
# The theoretical values are computed here independently of the code of SPOT, from       #
# the equations of the model (Supporting Information S1). Three checks:                  #
#   (a) power without drift and with perfect genotypes: exact calculation, by listing    #
#       all the possible allele counts at every sampling date;                           #
#   (b) neutral drift: expected genetic diversity left,                                  #
#       H_t / H_0 = product over generations of (1 - 1/(2 Ne_g))   (haploid: 1 - 1/N_g); #
#   (c) selection without drift: closed form for haploids,                               #
#       p_t = p0 (1+s)^t / (1 - p0 + p0 (1+s)^t), and the recursions of the model for    #
#       diploids and haplodiploids.                                                      #
#                                                                                        #
# Run from the repository root, after theory_spot.js:                                    #
#   python3 paper/validation/theory_check.py                                             #
# It writes theory_power.tsv, theory_drift.tsv, theory_selection.tsv and                 #
# theory_check.svg in paper/validation/, and prints a summary of each check.             #
##########################################################################################


#########################
# Load the dependencies #
#########################
import json, math, sys   # To read the JSON file; mathematical functions; Python search path
import numpy as np   # Arrays, used for the exact calculation of power
sys.path.insert(0, 'paper/validation')   # Folder that contains svgplot.py
from svgplot import Fig, Ax, COL   # Small plotting tools of this repository


######################################################
# Read the values given by SPOT and set fixed values #
######################################################
D = json.load(open('paper/validation/theory_spot.json'))   # Values given by SPOT (made by theory_spot.js)
Z95 = 1.6448536269514722   # Sampling-only threshold (95% quantile of the standard normal distribution)
TIE_TOL = 1e-9   # relative tolerance for ties at the threshold, as TIE_TOL in sim.js


######################################################################
# Function: allele frequencies under selection without drift         #
# Returns, for each generation g, the frequencies after selection:   #
# 'p' in females (or in all haploid individuals), 'pm' in males, and #
# 'pmf' the genotype frequencies of females (aa, Aa, AA)             #
######################################################################
def det(system, h, s, p0, G, smale=1.0):
    """deterministic recursion; element g = genotype/allele frequencies after selection in generation g (selection from g = 1)"""
    pf = pm = p0   # Frequency in females and in males at the start
    out = []   # Results, one element per generation
    for g in range(G + 1):   # Loop over the generations
        sg = s if g > 0 else 0.0   # Selection acts from generation 1
        sm = sg * smale   # Selection coefficient in males
        if system == 'haploid':
            ps = pf * (1 + sg) / (1 + pf * sg)   # Haploid: frequency after selection
            out.append({'p': ps, 'pm': ps, 'pmf': None})   # Keep the frequency of this generation (no genotypes for haploids)
            pf = pm = ps   # Frequency for the next generation
            continue   # Go to the next generation
        AA, Aa = pf * pm, pf * (1 - pm) + pm * (1 - pf)   # Genotype frequencies of the offspring (random union of a maternal and a paternal gene copy)
        aa = 1 - AA - Aa   # Frequency of the homozygotes without the allele
        w = AA * (1 + sg) + Aa * (1 + h * sg) + aa   # Mean fitness of females
        fAA, fAa, faa = AA * (1 + sg) / w, Aa * (1 + h * sg) / w, aa / w   # Genotype frequencies of females after selection
        pfs = fAA + fAa / 2   # Allele frequency in females after selection
        if system == 'haplodiploid':
            pms = pf * (1 + sm) / (pf * (1 + sm) + 1 - pf)   # Haplodiploid males: haploid, their gene copy comes from their mother
        else:
            wm = AA * (1 + sm) + Aa * (1 + h * sm) + aa   # Mean fitness of diploid males
            pms = (AA * (1 + sm) + Aa * (1 + h * sm) / 2) / wm   # Diploid males: allele frequency after selection
        out.append({'p': pfs, 'pm': pms, 'pmf': (faa, fAa, fAA)})   # Keep the frequencies of this generation
        pf, pm = pfs, pms   # Frequencies for the next generation
    return out   # Return the results


##################################################################
# Function: probability of each possible number of copies of the #
# allele among n_ind sampled individuals (females if diploid)    #
##################################################################
def count_pmf(system, state, n_ind):
    """distribution of the number of copies of the allele among n_ind sampled individuals (females if diploid)"""
    if system == 'haploid':
        p = state['p']   # Allele frequency in the population
        k = np.arange(n_ind + 1)   # Possible numbers of copies: 0 to n_ind
        return np.array([math.comb(n_ind, int(i)) * p ** i * (1 - p) ** (n_ind - i) for i in k])   # Haploid: binomial distribution
    g = np.array(state['pmf'])   # Diploid: probability that one female carries 0, 1 or 2 copies
    out = np.array([1.0])   # Start with no individual: 0 copies with probability 1
    for _ in range(n_ind):   # One step per sampled female
        out = np.convolve(out, g)   # Add one female at a time
    return out   # Return the results


######################################################################
# Function: score Z of the trend test for every possible combination #
# of allele counts at the sampling dates                             #
# nall: gene copies sampled at each date; t: generation of each date #
######################################################################
def z_grid(nall, t):
    """trend statistic for every combination of counts (x_1..x_K), alleles nall_i per point, times t_i"""
    K = len(nall)   # Number of sampling dates
    grids = np.meshgrid(*[np.arange(n + 1) for n in nall], indexing='ij')   # All the possible allele counts at each date
    N = sum(nall)   # Total number of gene copies sampled
    tb = sum(n * ti for n, ti in zip(nall, t)) / N   # Mean generation number of the sampled gene copies
    X = sum(grids)   # Total number of copies of the allele
    num = sum(x * (ti - tb) for x, ti in zip(grids, t))   # Numerator: covariance between allele copies and time
    den = sum(n * (ti - tb) ** 2 for n, ti in zip(nall, t))   # Sum of squared deviations of time
    pb = X / N   # Overall allele frequency in the samples
    with np.errstate(divide='ignore', invalid='ignore'):   # Do not warn about divisions by zero (handled just below)
        Z = num / np.sqrt(pb * (1 - pb) * den)   # Cochran-Armitage trend score
    Z[(pb <= 0) | (pb >= 1)] = 0.0   # Allele absent or fixed in all samples: Z = 0
    return Z   # Return Z for every combination


#########################################################################
# Function: joint probability of the allele counts of all dates         #
# (the dates are independent once the population frequencies are known) #
#########################################################################
def joint(pmfs):
    P = pmfs[0]   # Probabilities of the first date
    for q in pmfs[1:]:   # Add the other dates one by one
        P = np.multiply.outer(P, q)   # Product of the probabilities of each date
    return P   # Return the joint probabilities


###################################################################
# Function: exact threshold, power and false-positive rate of one #
# design of theory_spot.json, for the trend test and the          #
# first-versus-last test                                          #
###################################################################
def exact_power(r):
    system, h, s, p0, years, n = r['system'], r['h'], r['s'], r['p0'], r['years'], r['n']   # Settings of the design
    t = [y - years[0] for y in years]   # Generation of each sampling date (one generation per year)
    G = t[-1]   # Last generation
    sel = det(system, h, s, p0, G)   # Trajectory with selection
    neu = det(system, h, 0.0, p0, G)   # Trajectory without selection
    per = 1 if system == 'haploid' else 2   # Gene copies per sampled individual
    out = {}   # Results for the two tests
    for stat, idx in (('trend', list(range(len(t)))), ('two', [0, len(t) - 1])):   # Trend test (all dates), then first-versus-last test (first and last dates)
        nall = [n[i] * per for i in idx]   # Gene copies sampled at each date used
        tt = [t[i] for i in idx]   # Generation of each date used
        Z = z_grid(nall, tt)   # Z of every possible data set
        Ps = joint([count_pmf(system, sel[t[i]], n[i]) for i in idx])   # Probability of every possible data set with selection
        Pn = joint([count_pmf(system, neu[t[i]], n[i]) for i in idx])   # Probability of every possible data set without selection
        zs, pz = Z.ravel(), Pn.ravel()   # Z and its neutral probability, as flat lists
        o = np.argsort(zs, kind='stable')   # Order of Z from smallest to largest
        cdf = np.cumsum(pz[o])   # Cumulative distribution of Z without selection
        crit = zs[o][np.searchsorted(cdf, 0.95 - 1e-12)]   # smallest z with P(Z <= z) >= 0.95, as SPOT's empirical quantile
        # tie rule of sim.js (empiricalThreshold, exceed): values of Z equal to the threshold count with weight gamma.
        # Exact gamma at the exact quantile: P(Z > c) + gamma P(Z = c) = 0.05 exactly.
        eq = lambda c: np.abs(Z - c) <= TIE_TOL * max(1.0, abs(c))   # Function: is Z equal to c (a tie)?
        above = lambda P, c, g: float(P[(Z > c) & ~eq(c)].sum() + g * P[eq(c)].sum())   # Function: probability above c, ties counted with weight g
        A_, B_ = float(Pn[(Z > crit) & ~eq(crit)].sum()), float(Pn[eq(crit)].sum())   # Neutral probability strictly above the threshold, and of a tie
        gx = min(max((0.05 - A_) / B_, 0.0), 1 - 1e-12) if B_ > 0 else 0.0   # Weight of the ties that gives exactly 5% of false positives
        sc, gs = r[stat]['crit'], r[stat].get('gamma', 0.0)   # SPOT's own threshold and tie weight, from its neutral replicates
        out[stat] = {'crit': float(crit), 'gamma': gx, 'tieMass': B_, 'powerEmp': above(Ps, crit, gx), 'powerNaive': float(Ps[Z > Z95].sum()),
                     'fprNaive': float(Pn[Z > Z95].sum()), 'fprExact': above(Pn, crit, gx),
                     'powerAtSpot': above(Ps, sc, gs), 'sizeAtSpot': above(Pn, sc, gs)}
    return out   # Return the results


###############################################################
# (a) Power without drift: SPOT against the exact calculation #
###############################################################
rows, pts = [], []   # Lines of the results table, and points of the figure
for r in D['power']:   # Loop over the designs computed by SPOT
    ex = exact_power(r)   # Exact values for this design
    for stat in ('trend', 'two'):   # For both tests
        if r['design'] == '2 points' and stat == 'two':
            continue      # identical to the trend test with two points
        for kind, ekind in (('powerNaive', 'powerNaive'), ('powerEmp', 'powerAtSpot')):   # Sampling-only threshold, then the drift-aware threshold of SPOT
            a, b = ex[stat][ekind], r[stat][kind]   # Exact power, and power given by SPOT
            se = math.sqrt(max(a * (1 - a), 1e-12) / r['reps'])   # Standard error expected for the estimate of SPOT
            z = (b - a) / se if a * (1 - a) > 0 else 0.0   # Difference between SPOT and the exact value, in standard errors
            rows.append([r['system'], r['h'], r['p0'], r['s'], r['design'], stat, kind, f'{a:.4f}', f'{b:.4f}', f'{z:.2f}',
                         f'{ex[stat]["crit"]:.3f}', f'{r[stat]["crit"]:.3f}', f'{ex[stat]["fprExact"]:.4f}', f'{ex[stat]["sizeAtSpot"]:.4f}',
                         f'{ex[stat]["powerEmp"]:.4f}'])
            pts.append((a, b, kind, stat, z, ex[stat]['crit'], r[stat]['crit'], r['system'], ex[stat]['sizeAtSpot'], ex[stat]['fprExact'], ex[stat]['powerEmp']))
with open('paper/validation/theory_power.tsv', 'w') as f:   # Write the table of the power check
    f.write('\t'.join(['system', 'h', 'p0', 's', 'design', 'statistic', 'threshold', 'power_exact_at_spot_threshold', 'power_spot', 'z',
                       'crit_exact', 'crit_spot', 'size_at_exact_crit', 'size_at_spot_crit', 'power_at_exact_crit']) + '\n')
    for r in rows:   # One line per comparison
        f.write('\t'.join(map(str, r)).replace('powerNaive', 'naive z(0.95)').replace('powerEmp', 'neutral 95% quantile') + '\n')
zs = np.array([p[4] for p in pts])   # Differences in standard errors, all comparisons
d = np.array([p[1] - p[0] for p in pts])   # Differences SPOT minus exact, all comparisons
crit_same = sum(abs(p[5] - p[6]) < 1e-6 for p in pts)   # Number of cases where the threshold of SPOT equals the exact one
print(f'(a) power: {len(pts)} comparisons; mean SPOT - exact {d.mean() * 100:+.2f} points; max |diff| {abs(d).max() * 100:.1f} points;'
      f' |z| > 2: {(abs(zs) > 2).sum()}, > 3: {(abs(zs) > 3).sum()}; RMS z {math.sqrt((zs ** 2).mean()):.2f};'
      f' threshold identical to exact in {crit_same} of {len(pts)}')
for p in pts:   # Print the comparisons that differ by more than 3 standard errors
    if abs(p[4]) > 3:
        print('   large z:', p)
for kind in ('powerNaive', 'powerEmp'):   # Summary for each threshold
    zz = np.array([p[4] for p in pts if p[2] == kind])   # Differences in standard errors for this threshold
    dd = np.array([p[1] - p[0] for p in pts if p[2] == kind])   # Differences SPOT minus exact for this threshold
    print(f'    {kind}: n {len(zz)}, mean diff {dd.mean()*100:+.2f}, max |diff| {abs(dd).max()*100:.1f}, |z|>2 {(abs(zz)>2).sum()}, RMS z {math.sqrt((zz**2).mean()):.2f}')
em = [p for p in pts if p[2] == 'powerEmp']   # Comparisons with the drift-aware threshold
sz = np.array([p[8] for p in em])   # Exact false-positive rate of the thresholds of SPOT
sx = np.array([p[9] for p in em])   # Exact false-positive rate of the exact thresholds
print(f'    size of SPOT thresholds (exact): mean {sz.mean():.4f}, range {sz.min():.4f}-{sz.max():.4f}; size at exact 95% quantile: mean {sx.mean():.4f}, range {sx.min():.4f}-{sx.max():.4f}')
dpow = np.array([p[1] - p[10] for p in em])   # Power of SPOT minus exact power at the exact threshold
print(f'    SPOT power vs exact power at the exact quantile: mean {dpow.mean()*100:+.2f}, max |diff| {abs(dpow).max()*100:.1f} points')


##################################################################
# (b) Neutral drift: genetic diversity left, SPOT against theory #
##################################################################
def ne(system, Nf, Nm):   # Function: effective population size from the breeding numbers
    if system == 'haploid':
        return Nf   # Haploid: N
    if system == 'haplodiploid':
        return 9 * Nf * Nm / (4 * Nm + 2 * Nf)   # Haplodiploid
    return 4 * Nf * Nm / (Nf + Nm)   # Diploid


#########################################################################
# Compare the simulated diversity with theory, generation by generation #
#########################################################################
drift_rows, drift_z = [], []   # Results for the table and the figure
for r in D['drift']:   # Loop over the drift scenarios
    o, G = r['params'], r['G']   # Settings of the scenario, and its last generation
    sysm, Nf0, Nm0 = o['system'], o['Nf'], o['Nm']   # Inheritance system and breeding numbers at the start
    gens = [0, 20, 40, 90]          # sampling generations of the default design (1935, 1955, 1975, 2025)
    H, theo = 1.0, [1.0]   # Diversity left according to theory, starting at 1
    for g in range(G):   # Loop over the generations
        if o.get('neByBlock'):   # Ne set separately for each interval between sampling dates
            j = 0   # Find the interval that contains generation g
            while j < len(gens) - 2 and g >= gens[j + 1]:
                j += 1
            f = o['neBlocks'][j] / ne(sysm, Nf0, Nm0)   # Ne of this interval relative to the Ne given by Nf and Nm
        else:
            f = o.get('decline', 1) ** (g / G)   # Exponential change of the population size
        Nf, Nm = max(2, round(Nf0 * f)), max(2, round(Nm0 * f))   # Breeding numbers of this generation
        e = ne(sysm, Nf, Nm)   # Ne of this generation
        H *= (1 - 1 / e) if sysm == 'haploid' else (1 - 1 / (2 * e))   # Diversity lost by drift in one generation
        theo.append(H)   # Keep the value of this generation
    sim, se = r['het'], r['se']   # Simulated diversity left and its standard error (from SPOT)
    rel = [sim[g] / theo[g] - 1 for g in range(1, G + 1)]   # Relative difference between simulation and theory
    drift_rows.append((r['label'], theo, sim, r['spotHetRetained']))   # Keep the results of this scenario
    zmax = max(abs(sim[g] - theo[g]) / se[g] for g in range(1, G + 1))   # Largest difference over the generations, in standard errors
    drift_z.append((r['label'], (sim[-1] - theo[-1]) / se[-1]))   # Difference at the last generation, in standard errors
    print(f'(b) {r["label"]:45s} H90/H0 sim {sim[-1]:.4f} ± {se[-1]:.4f} theory {theo[-1]:.4f} (SPOT tile {r["spotHetRetained"]:.4f});'
          f' z at generation 90 {(sim[-1] - theo[-1]) / se[-1]:+.2f}; max relative deviation {max(abs(x) for x in rel) * 100:.2f}%, max |z| over generations {zmax:.2f}')
with open('paper/validation/theory_drift.tsv', 'w') as f:   # Write the table of the drift check
    f.write('scenario\tgeneration\tH_ratio_simulated\tH_ratio_theory\n')
    for lab, th, sm, _ in drift_rows:   # One block of lines per scenario
        for g in range(len(th)):   # One line per generation
            f.write(f'{lab}\t{g}\t{sm[g]:.5f}\t{th[g]:.5f}\n')


####################################################################
# (c) Selection without drift: trajectories of SPOT against theory #
####################################################################
maxdiff, sel_rows = 0.0, []   # Largest difference over all trajectories, and results per trajectory
for r in D['selection']:   # Loop over the trajectories computed by SPOT
    G = len(r['f']) - 1   # Last generation
    if r['system'] == 'haploid':
        th = [r['p0'] * (1 + r['s']) ** g / (1 - r['p0'] + r['p0'] * (1 + r['s']) ** g) for g in range(G + 1)]   # Haploid: closed form
        thm = th   # Haploid: no separate value for males
    else:
        dd = det(r['system'], r['h'], r['s'], r['p0'], G)   # Diploid and haplodiploid: recursion of the model
        th, thm = [x['p'] for x in dd], [x['pm'] for x in dd]
    md = max(max(abs(a - b) for a, b in zip(r['f'], th)), max(abs(a - b) for a, b in zip(r['m'], thm)))   # Largest difference between SPOT and theory along the trajectory
    maxdiff = max(maxdiff, md)   # Keep the largest difference so far
    sel_rows.append((r, th, md))   # Keep the results of this trajectory
print(f'(c) selection: {len(sel_rows)} trajectories (300 generations); largest |SPOT - theory| = {maxdiff:.2e}')
sm = D['stochasticMean']   # Mean of 500 trajectories with drift in a very large population (from SPOT)
dd = det('diploid', 0.5, sm['s'], sm['p0'], 300)   # The same trajectory without drift, from theory
print(f'    stochastic mean of {sm["reps"]} replicates (N = 1e6) vs deterministic: max |diff| = '
      f'{max(abs(a - b["p"]) for a, b in zip(sm["f"], dd)):.4f}')
with open('paper/validation/theory_selection.tsv', 'w') as f:   # Write the table of the selection check
    f.write('system\th\ts\tp0\tmax_abs_difference_spot_vs_theory\n')
    for r, _, md in sel_rows:   # One line per trajectory
        f.write(f'{r["system"]}\t{r["h"]}\t{r["s"]}\t{r["p0"]}\t{md:.2e}\n')


#################################################
# Figure of the three checks (theory_check.svg) #
#################################################
F = Fig(1180, 440)   # Empty figure (width and height in pixels)
# (a) power scatter
A = Ax(F, 70, 50, 300, 300, (0, 1), (0, 1), title='Power without drift: SPOT vs exact', label='a')   # Panel a: power of SPOT against exact power
A.frame([0, .2, .4, .6, .8, 1], [0, .2, .4, .6, .8, 1], 'exact power', 'SPOT power')   # Ticks and axis titles
xs = [i / 100 for i in range(101)]   # Values from 0 to 1, to draw the grey band
se = [2 * math.sqrt(x * (1 - x) / 5000) for x in xs]   # Two standard errors of a power estimated from 5,000 simulations
A.band(xs, [x - e for x, e in zip(xs, se)], [x + e for x, e in zip(xs, se)], COL['grey'], .25)   # Grey band: +/- 2 standard errors around the 1:1 line
A.line([0, 1], [0, 1], COL['muted'], 1, '4 4')   # Dashed 1:1 line (perfect agreement)
for p in pts:   # Print the comparisons that differ by more than 3 standard errors
    A.points([p[0]], [p[1]], COL['blue'] if p[2] == 'powerEmp' else COL['red'], 'o' if p[3] == 'trend' else 's', 3.2)   # Blue: drift-aware threshold, red: sampling-only; circle: trend test, square: first-versus-last
A.legend([('drift-aware threshold (SPOT\'s own)', COL['blue'], 'o'), ('sampling-only threshold z₀.₉₅', COL['red'], 'o'),
          ('trend test, all points', COL['muted'], 'o'), ('first vs last point', COL['muted'], 's')], A.x + 8, A.y + 16)
# (b) heterozygosity
B = Ax(F, 470, 50, 300, 300, (0, 90), (0, 1), title='Neutral drift: heterozygosity retained', label='b')   # Panel b: genetic diversity left by drift
B.frame([0, 20, 40, 60, 80], [0, .2, .4, .6, .8, 1], 'generation', 'H(t) / H(0)')   # Ticks and axis titles
pal = [COL['blue'], COL['green'], COL['amber'], COL['purple'], COL['red']]   # Colours used in turn for the curves
show = ['haploid, N = 100', 'diploid, Nf = Nm = 50', 'haplodiploid, Nf = Nm = 100', 'diploid, Nf = Nm = 500, decline to 10%',
        'diploid, Ne per interval 1,000 / 50 / 1,000']
leg = []   # Entries of the legend
for k, lab in enumerate(show):   # One curve per scenario shown
    th, sim = next((t, s) for l, t, s, _ in drift_rows if l == lab)   # Theory and simulation of this scenario
    B.line(list(range(len(th))), th, pal[k], 1.6)   # Line: theory
    g = list(range(0, len(sim), 6))   # One simulated point every 6 generations
    B.points(g, [sim[i] for i in g], pal[k], 'o', 3, True)   # Hollow circles: SPOT
    leg.append((lab.replace('Ne per interval', 'Ne'), pal[k], 'line'))   # Legend entry of this scenario
B.legend(leg + [('lines: theory, circles: SPOT', COL['muted'], 'oh')], B.x + 8, B.y + B.h - 100)   # Legend, at the bottom left of the panel
# (c) haploid selection
C = Ax(F, 870, 50, 280, 300, (0, 300), (0, 1), title='Selection without drift (haploid)', label='c')   # Panel c: haploid trajectories under selection without drift
C.frame([0, 50, 100, 150, 200, 250, 300], [0, .2, .4, .6, .8, 1], 'generation', 'allele frequency')   # Ticks and axis titles
for k, s in enumerate([0.01, 0.05, 0.2]):   # One colour per selection coefficient
    for p0 in (0.01, 0.2):   # Two starting frequencies
        r, th, _ = next(x for x in sel_rows if x[0]['system'] == 'haploid' and x[0]['s'] == s and x[0]['p0'] == p0)   # SPOT trajectory and theory for this case
        C.line(list(range(0, 301, 2)), th[::2], pal[k], 1.6)   # Line: closed form
        g = list(range(0, 301, 10))   # One SPOT point every 10 generations
        C.points(g, [r['f'][i] for i in g], pal[k], 'o', 2.6, True)   # Hollow circles: SPOT
C.legend([('s = 0.01', pal[0], 'line'), ('s = 0.05', pal[1], 'line'), ('s = 0.2', pal[2], 'line')], C.X(185), C.Y(0.42))   # Legend
F.text(70, 420, f'a: {len(D["power"])} designs × 2 thresholds; shaded: ±2 Monte Carlo SE (5,000 replicates). b: 100,000 neutral trajectories each, p₀ = 0.3. '   # Notes under the figure
       'c: p₀ = 0.01 and 0.2; lines: closed form, circles: SPOT.', 11, color=COL['muted'])
print('figure bytes', F.save('paper/validation/theory_check.svg'))   # Save the figure
