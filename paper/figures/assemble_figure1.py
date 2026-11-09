"""Assembles Figure 1 from the five panels saved by make_figure1.js, in two columns of equal width:
a, result tiles; b, distribution of Z; c, power curve; d, sample-size planner, with e, allele trajectories, below d.
Every panel is scaled to the same
width, so text has about the same size everywhere. Panel letters top left."""
import sys, os
from PIL import Image, ImageDraw, ImageFont
D = sys.argv[1]
P = {k: Image.open(os.path.join(D, 'Figure1%s.png' % k)).convert('RGB') for k in 'abcde'}
Wc, G, TOP = 1600, 70, 90
def fit(im): return im.resize((Wc, round(im.height * Wc / im.width)), Image.LANCZOS)
P = {k: fit(v) for k, v in P.items()}
h1 = max(P['a'].height, P['b'].height); h2 = max(P['c'].height, P['d'].height + TOP + P['e'].height)
W = 2 * Wc + G; H = TOP + h1 + G + TOP + h2 + 20
F = Image.new('RGB', (W, H), 'white')
pos = {'a': (0, TOP), 'b': (Wc + G, TOP), 'c': (0, TOP + h1 + G + TOP), 'd': (Wc + G, TOP + h1 + G + TOP),
       'e': (Wc + G, TOP + h1 + G + TOP + P['d'].height + TOP)}
dr = ImageDraw.Draw(F)
try: font = ImageFont.truetype('/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf', 64)
except OSError: font = ImageFont.load_default()
for k, (x, y) in pos.items():
    F.paste(P[k], (x, y)); dr.text((x, y - 80), k, fill='black', font=font)
F.save(os.path.join(D, 'Figure1.png'), dpi=(300, 300))
print('Figure1.png', F.size)
