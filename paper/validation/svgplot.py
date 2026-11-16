####################################################################################
# Small plotting tools used by the figure scripts of this repository               #
#                                                                                  #
# The figures are written directly as SVG text, with the Python standard library   #
# only. They stay small enough to be viewed on GitHub, and they can be reproduced  #
# without installing a plotting package.                                           #
#                                                                                  #
# How a figure script uses this file:                                              #
#   F = Fig(width, height)                     an empty figure                     #
#   a = Ax(F, x, y, w, h, xlim, ylim, ...)     a panel (axes) placed in the figure #
#   a.frame(...); a.line(...); a.points(...)   draw in the panel                   #
#   F.save('file.svg')                         write the SVG file                  #
# Positions and sizes are in pixels; the origin is the top left corner.            #
####################################################################################


#########################
# Load the dependencies #
#########################
import math                                   # For logarithms (log scales, powers of ten)


##########################
# Colours and text sizes #
##########################
COL = {'blue': '#2B59D9', 'green': '#0E9F6E', 'red': '#D9363E', 'amber': '#C27A12', 'purple': '#8B4FD6',
       'grey': '#9AA3B8', 'ink': '#0F1730', 'muted': '#5A6480', 'line': '#E1E5EF'}   # Colours shared by all figures (ink: text; muted: ticks; line: grid)
# text sizes of the axes; a figure script can change them before drawing (see fpr_figure.py)
STYLE = {'tick': 11, 'lab': 12, 'title': 13, 'title_centred': False}   # Sizes of tick numbers, axis titles and panel titles; panel title centred or not


####################################################
# Function: make a text safe to put in an SVG file #
####################################################
def esc(t):
    return str(t).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')   # Replace the three characters that have a special meaning in SVG


#################################################################################
# Function: write a number as a tick label                                      #
# 10,000 and above: power of ten (10⁴, 2×10⁵); otherwise thousands separator or #
# plain number                                                                  #
#################################################################################
def fmt(v):
    if v == 0:
        return '0'
    a = abs(v)
    if a >= 10000:                            # Large number: written as m×10^e
        e = int(math.floor(math.log10(a) + 1e-9))   # Exponent
        m = v / 10 ** e                             # Number in front of the power of ten
        return ('10' if abs(m - 1) < 1e-9 else f'{m:g}×10') + f'<tspan dy="-5" font-size="{round(STYLE["tick"] * 8 / 11)}">{e}</tspan>'   # The exponent is raised and smaller
    if a >= 1:
        return f'{v:,.0f}' if abs(v - round(v)) < 1e-9 else f'{v:g}'   # Whole number: with a thousands separator
    return f'{v:g}'


#############################################################################
# Class Fig: a figure                                                       #
# It keeps the list of SVG elements (self.o) and writes the file at the end #
#############################################################################
class Fig:
    def __init__(self, W, H, title=None):
        self.W, self.H, self.o = W, H, [f'<rect width="{W}" height="{H}" fill="#fff"/>']   # Width, height, and the list of elements (starts with a white background)
        self.defs, self.defo = {}, []         # Marker shapes already defined (see Ax.marker)
        if title:                             # Optional title at the top of the figure
            self.o.append(f'<text x="{W / 2}" y="22" text-anchor="middle" font-size="15" font-weight="600" fill="{COL["ink"]}">{esc(title)}</text>')

    # Add a text at position (x, y): size, alignment ('start', 'middle' or 'end'), colour, weight, rotation in degrees
    def text(self, x, y, t, size=12, anchor='start', color=None, weight=None, rot=None):
        tr = f' transform="rotate({rot} {x:.1f} {y:.1f})"' if rot else ''   # Rotation around the position of the text
        w = f' font-weight="{weight}"' if weight else ''
        self.o.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" text-anchor="{anchor}" fill="{color or COL["ink"]}"{w}{tr}>{t}</text>')

    # Write the SVG file and return its size in characters
    def save(self, path):
        s = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.W}" height="{self.H}" viewBox="0 0 {self.W} {self.H}" '
             f'font-family="Helvetica,Arial,sans-serif">\n<defs>' + ''.join(self.defo) + '</defs>\n' + '\n'.join(self.o) + '\n</svg>\n')
        open(path, 'w').write(s)
        return len(s)


##################################################################################
# Class Ax: a panel (a pair of axes) placed in a figure                          #
# x, y: top left corner; w, h: width and height; xlim, ylim: limits of the axes; #
# xlog, ylog: log scales; title: panel title; label: panel letter                #
##################################################################################
class Ax:
    def __init__(self, fig, x, y, w, h, xlim, ylim, xlog=False, ylog=False, title=None, label=None):
        self.f, self.x, self.y, self.w, self.h = fig, x, y, w, h
        self.xlim, self.ylim, self.xlog, self.ylog = xlim, ylim, xlog, ylog
        if title:                             # Panel title, above the panel
            if STYLE['title_centred']:
                fig.text(x + w / 2, y - 10, esc(title), STYLE['title'], 'middle', weight='600')   # Centred above the panel
            else:
                fig.text(x, y - 10, esc(title), STYLE['title'], weight='600')                     # Aligned with the left of the panel
        if label:                             # Panel letter, at the top left
            fig.text(x - 52, y - 10, label, 15, weight='700')

    # Position of the value v between the two limits of an axis, from 0 to 1 (log scale if asked)
    def _t(self, v, lim, log):
        a, b = lim
        if log:
            v, a, b = math.log10(v), math.log10(a), math.log10(b)
        return (v - a) / (b - a)

    # Horizontal position (pixels) of the x value v
    def X(self, v):
        return self.x + self._t(v, self.xlim, self.xlog) * self.w

    # Vertical position (pixels) of the y value v (the y axis of an SVG points down, hence the minus)
    def Y(self, v):
        return self.y + self.h - self._t(v, self.ylim, self.ylog) * self.h

    # Draw the frame of the panel: grid, ticks, tick numbers and axis titles
    # xt, yt: tick values; xlab, ylab: axis titles; xfmt, yfmt: functions that write a tick number
    def frame(self, xt, yt, xlab='', ylab='', xfmt=fmt, yfmt=fmt, grid=True):
        o, c = self.f.o, COL
        ts, ls = STYLE['tick'], STYLE['lab']   # tick and axis-label sizes; the labels move out as the text grows
        for v in yt:                          # y axis: one grid line and one number per tick
            yy = self.Y(v)
            if grid:
                o.append(f'<line x1="{self.x}" x2="{self.x + self.w}" y1="{yy:.1f}" y2="{yy:.1f}" stroke="{c["line"]}"/>')
            o.append(f'<text x="{self.x - 6}" y="{yy + ts * 4 / 11:.1f}" font-size="{ts}" text-anchor="end" fill="{c["muted"]}">{yfmt(v)}</text>')
        for v in xt:                          # x axis: one small tick and one number per tick
            xx = self.X(v)
            o.append(f'<line x1="{xx:.1f}" x2="{xx:.1f}" y1="{self.y + self.h}" y2="{self.y + self.h + 4}" stroke="{c["muted"]}"/>')
            o.append(f'<text x="{xx:.1f}" y="{self.y + self.h + 6 + ts}" font-size="{ts}" text-anchor="middle" fill="{c["muted"]}">{xfmt(v)}</text>')
        o.append(f'<rect x="{self.x}" y="{self.y}" width="{self.w}" height="{self.h}" fill="none" stroke="{c["ink"]}" stroke-width=".8"/>')   # Rectangle around the panel
        if xlab:                              # Title of the x axis, under the tick numbers
            self.f.text(self.x + self.w / 2, self.y + self.h + 36 + 2 * (ts - 11) + (ls - 12), xlab, ls, 'middle', c['ink'])
        if ylab:                              # Title of the y axis, written vertically left of the tick numbers
            self.f.text(self.x - 44 - 4 * (ts - 11), self.y + self.h / 2, ylab, ls, 'middle', c['ink'], rot=-90)

    # Draw a line through the points (xs, ys); dash: pattern of a dashed line, for example '5 4'
    def line(self, xs, ys, color, width=1.8, dash=None, opacity=1):
        d = 'M' + 'L'.join(f'{self.X(a):.1f} {self.Y(b):.1f}' for a, b in zip(xs, ys))   # Path: move to the first point, then lines to the next ones
        da = f' stroke-dasharray="{dash}"' if dash else ''
        self.f.o.append(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{width}"{da} stroke-opacity="{opacity}"/>')

    # Draw a horizontal dashed line at the y value v, across the whole panel
    def hline(self, v, color, dash='5 4', width=1.2):
        self.line([self.xlim[0], self.xlim[1]], [v, v], color, width, dash)

    # Draw a shaded band between the curves lo and hi
    def band(self, xs, lo, hi, color, opacity=.15):
        pts = [(self.X(a), self.Y(b)) for a, b in zip(xs, hi)] + [(self.X(a), self.Y(b)) for a, b in zip(reversed(xs), reversed(lo))]   # Along the upper curve, then back along the lower one
        self.f.o.append('<path d="M' + 'L'.join(f'{a:.1f} {b:.1f}' for a, b in pts) + f'Z" fill="{color}" fill-opacity="{opacity}"/>')

    # Return the SVG text of one marker at the pixel position (px, py)
    # shape: 'o' circle, 's' square, '^' triangle, anything else diamond; hollow: white inside
    def marker(self, px, py, color, shape='o', size=4, hollow=False):
        """markers are defined once in <defs> and placed with <use>, which keeps files small"""
        key = (color, shape, size, hollow)    # A marker shape is defined once for each combination
        f = self.f
        if key not in f.defs:                 # First use of this combination: define the shape
            mid = f'm{len(f.defs)}'           # Name of the shape (m0, m1, ...)
            fill, st = ('#fff', color) if hollow else (color, '#fff')   # Hollow: white inside, coloured outline; filled: the reverse
            sw = 1.5 if hollow else .6        # Width of the outline
            s = size
            if shape == 'o':
                g = f'<circle r="{s}"'
            elif shape == 's':
                g = f'<rect x="{-s}" y="{-s}" width="{2 * s}" height="{2 * s}"'
            elif shape == '^':
                g = f'<path d="M0 {-s - 1}L{s + .5} {s - .5}H{-s - .5}Z"'
            else:
                g = f'<path d="M0 {-s - 1}L{s + 1} 0L0 {s + 1}L{-s - 1} 0Z"'
            f.defs[key] = mid
            f.defo.append(f'<g id="{mid}">{g} fill="{fill}" stroke="{st}" stroke-width="{sw}"/></g>')
        return f'<use href="#{f.defs[key]}" x="{px:.1f}" y="{py:.1f}"/>'   # Place a copy of the shape

    # Draw one marker at each point (xs, ys)
    def points(self, xs, ys, color, shape='o', size=4, hollow=False):
        for a, b in zip(xs, ys):
            self.f.o.append(self.marker(self.X(a), self.Y(b), color, shape, size, hollow))

    # Draw a vertical bar from lo to hi at each x value (error bars)
    def errbar(self, xs, lo, hi, color, width=1.2):
        for a, l, h in zip(xs, lo, hi):
            xx = self.X(a)
            self.f.o.append(f'<line x1="{xx:.1f}" x2="{xx:.1f}" y1="{self.Y(l):.1f}" y2="{self.Y(h):.1f}" stroke="{color}" stroke-width="{width}"/>')

    # Draw a legend: one line per item, starting at the pixel position (x, y); dy: distance between lines
    def legend(self, items, x=None, y=None, dy=17):
        """items: (label, color, kind) with kind 'line', 'dash', or a marker shape ('o', 's', '^', 'd'); 'o-' = line + marker"""
        x = self.x + 10 if x is None else x   # Default position: top left of the panel
        y = self.y + 14 if y is None else y
        for i, (lab, color, kind) in enumerate(items):
            yy = y + i * dy                   # Vertical position of this item
            if kind in ('line', 'dash', 'o-', 's-', '^-'):                # Item drawn as a line (with a marker if the kind ends with '-')
                da = ' stroke-dasharray="5 3"' if kind == 'dash' else ''
                self.f.o.append(f'<line x1="{x}" x2="{x + 20}" y1="{yy - 4}" y2="{yy - 4}" stroke="{color}" stroke-width="2"{da}/>')
                if kind.endswith('-'):
                    self.f.o.append(self.marker(x + 10, yy - 4, color, kind[0], 3.5))
            else:                                                         # Item drawn as a marker only (hollow if the kind ends with 'h')
                self.f.o.append(self.marker(x + 10, yy - 4, color, kind.rstrip('h'), 4, kind.endswith('h')))
            self.f.text(x + 27, yy, esc(lab), 11, color=COL['ink'])       # Text of the item
