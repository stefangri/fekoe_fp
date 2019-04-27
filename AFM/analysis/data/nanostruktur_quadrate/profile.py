import numpy as np
import uncertainties.unumpy as unp
import matplotlib as mlp
mlp.use('pgf')
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.image as mpimg

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = '[lmroman10-regular]:+tlig;'
rcParams['text.usetex'] = True
rcParams['pgf.texsystem'] = 'lualatex'
rcParams['font.size'] = 10
rcParams['mathtext.fontset'] = 'custom'
rcParams['pgf.preamble'] = r'\usepackage[locale=DE]{siunitx}'#\usepackage{tgpagella}
rcParams['text.latex.preamble'] = r'\usepackage[math-style=ISO,bold-style=ISO,sans-style=italic,nabla=upright,partial=upright,]{unicode-math}'
#rcParams['text.latex.preamble'] = r'\usepackage[locale=DE]{siunitx}'
#rcParams['text.latex.preamble'] = r'\DeclareSIUnit{\counts}{Counts}'
rcParams['axes.formatter.use_mathtext'] = True
rcParams['legend.fontsize'] = 10
rcParams['figure.figsize'] =  5.906, 5.906
rcParams['savefig.dpi'] = 300

x1, y1, x2, y2 = np.genfromtxt('quadrate_profile.txt', unpack = True)
x1 *= 1e6
y1 *= 1e6
x2 *= 1e6
y2 *= 1e6

fig = plt.figure()

#image = mpimg.imread('quadrate_3d.png')
#ax0 = fig.add_subplot(222)
#
#ax0.imshow(image[:, 1000:], interpolation='nearest')
#ax0.axis('off')

image = mpimg.imread('quadrate_250_250_100pps.png')
ax3 = fig.add_subplot(211)
ax3.imshow(image)
ax3.set_yticks([])
ax3.set_xticks([])
ax3.text(0.01 , 0.96, r'(a)', horizontalalignment='left', verticalalignment='top', transform = ax3.transAxes, color = 'white')

ax1 = fig.add_subplot(223)

ax1.plot(x1, y1, label = 'Profil horizontal')
ax1.set_xlim(x1[0], x1[-8])
ax1.set_xlabel(r'Position $x / \si{\micro\meter}$ ')
ax1.set_ylabel('Höhe $y / \si{\micro\meter}$')
ax1.legend(loc = 'lower right')
ax1.text(0.01 , 0.96, r'(b)', horizontalalignment='left', verticalalignment='top', transform = ax1.transAxes)


ax2 = fig.add_subplot(224, sharey = ax1)

ax2.plot(x2, y2, label = 'Profil vertikal')
ax2.set_xlim(x2[0], x2[-8])
ax2.set_xlabel(r'Position $x / \si{\micro\meter}$ ')
ax2.legend(loc = 'lower right')
ax2.text(0.01 , 0.96, r'(c)', horizontalalignment='left', verticalalignment='top', transform = ax2.transAxes)
#ax2.set_ylabel('Höhe $y / \si{\micro\meter}$')


fig.tight_layout()
fig.savefig('quadrate_profile.pdf', bbox_inches = 'tight', pad_inches = 0)
