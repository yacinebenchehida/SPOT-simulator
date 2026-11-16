#################################################################################
# Assemble Figure 1 from the five panels saved by make_figure1.js               #
#                                                                               #
# Layout: two columns of equal width.                                           #
#   left  : a (result tiles), then c (power curve)                              #
#   right : b (distribution of Z), then d (sample-size planner) with e (allele  #
#           trajectories) below d                                               #
# Every panel is scaled to the same width, so that the text has about the same  #
# size everywhere. The panel letters are written at the top left of each panel. #
#                                                                               #
# Called by make_figure1.js. To run it alone, from the repository root:         #
#   python3 paper/figures/assemble_figure1.py paper/figures                     #
#################################################################################


#########################
# Load the dependencies #
#########################
import sys                                    # To read the argument of the command line
import os                                     # To build file names
from PIL import Image, ImageDraw, ImageFont   # Pillow: to open, resize, paste and write on images


########################
# Read the five panels #
########################
D = sys.argv[1]                               # Folder that contains Figure1a.png ... Figure1e.png
P = {k: Image.open(os.path.join(D, 'Figure1%s.png' % k)).convert('RGB') for k in 'abcde'}   # The five panels, by letter


#######################################
# Scale every panel to the same width #
#######################################
Wc, G, TOP = 1600, 70, 90                     # Width of a column, gap between columns and rows, space above a panel for its letter (pixels)

def fit(im):
    """Return the image resized to the width of a column, keeping its proportions."""
    return im.resize((Wc, round(im.height * Wc / im.width)), Image.LANCZOS)

P = {k: fit(v) for k, v in P.items()}


######################
# Size of the figure #
######################
h1 = max(P['a'].height, P['b'].height)                               # Height of the first row (a and b)
h2 = max(P['c'].height, P['d'].height + TOP + P['e'].height)         # Height of the second row (c on the left; d and e on the right)
W = 2 * Wc + G                                                       # Width of the figure
H = TOP + h1 + G + TOP + h2 + 20                                     # Height of the figure
F = Image.new('RGB', (W, H), 'white')                                # Empty white figure


############################################
# Position of each panel (top left corner) #
############################################
pos = {'a': (0, TOP),
       'b': (Wc + G, TOP),
       'c': (0, TOP + h1 + G + TOP),
       'd': (Wc + G, TOP + h1 + G + TOP),
       'e': (Wc + G, TOP + h1 + G + TOP + P['d'].height + TOP)}      # e is below d


############################################
# Paste the panels and write their letters #
############################################
dr = ImageDraw.Draw(F)                        # Tool to write on the figure
try:
    font = ImageFont.truetype('/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf', 64)   # Bold font for the letters
except OSError:
    font = ImageFont.load_default()           # Default font if the bold font is not installed
for k, (x, y) in pos.items():
    F.paste(P[k], (x, y))                     # Paste the panel
    dr.text((x, y - 80), k, fill='black', font=font)   # Write its letter above its top left corner


###################
# Save the figure #
###################
F.save(os.path.join(D, 'Figure1.png'), dpi=(300, 300))
print('Figure1.png', F.size)                  # Print the size of the figure in pixels
