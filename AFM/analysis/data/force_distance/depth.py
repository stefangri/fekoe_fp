import numpy as np
import uncertainties.unumpy as unp
import matplotlib as mlp
mlp.use('pgf')
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.image as mpimg
from scipy.optimize import curve_fit
from uncertainties import correlated_values
from uncertainties import ufloat

from pint import UnitRegistry
u = UnitRegistry()
Q_ = u.Quantity

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
rcParams['figure.figsize'] =  5.906, 2.5
rcParams['savefig.dpi'] = 300


xt, yt = np.genfromtxt('teflon_repulsive.txt', unpack = True)
xes, yes = np.genfromtxt('edelstahl_repulsive.txt', unpack = True)
xdlc, ydlc = np.genfromtxt('dlc_repulsive.txt', unpack = True)
xt -= xt[0]
yt -= yt[0]
xes -= xes[0]
yes -= yes[0]
xdlc -= xdlc[0]
ydlc -= ydlc[0]

fig = plt.figure()

ax1 = fig.add_subplot(131)

ax1.plot(xt, yt , label = 'Teflon')
ax1.plot(xes[xes <= xt[-1]], yes[xes <= xt[-1]], label = 'Edelstahl', color = 'g')
ax1.set_ylabel(r'Differenzspannung / \si{\volt}')

depth = xt[-1] - (xes[yes <= yt[-1]])[-1]
ax1.axvline(x = xt[-1], color = 'grey', linestyle = '--')
ax1.axvline(x = xt[-1] - depth, color = 'grey', linestyle = '--')
poisson = 0.46 #@23 Grad
k = Q_(0.2, 'newton / meter')
alpha = 10 / 180 * np.pi

depth = Q_(depth, 'micrometer')
print(depth)
z = Q_(xt[-1], 'micrometer')
E = k * (z - depth) * np.pi * (1 - poisson**2) / (2 * np.tan(alpha) * depth**2)
print(E.to('kPa'))




ax2 = fig.add_subplot(132)
ax2.plot(xt[xt <= xdlc[-1]], yt[xt <= xdlc[-1]], label = 'Teflon')
ax2.plot(xdlc, ydlc, label = 'DLC', color = 'r')

depth = (xt[xt <= xdlc[-1]])[-1] - (xdlc[ydlc <= (yt[xt <= xdlc[-1]])[-1]])[-1]
ax2.axvline(x = (xt[xt <= xdlc[-1]])[-1], color = 'grey', linestyle = '--')
ax2.axvline(x = (xt[xt <= xdlc[-1]])[-1] - depth, color = 'grey', linestyle = '--')
ax2.set_xlabel('Verschiebung $\Delta x / \si{\micro\meter}$')
print(depth)

depth = Q_(depth, 'micrometer')
z = Q_((xt[xt <= xdlc[-1]])[-1], 'micrometer')
E = k * (z - depth) * np.pi * (1 - poisson**2) / (2 * np.tan(alpha) * depth**2)
print(E.to('kPa'))

ax3 = fig.add_subplot(133)
ax3.plot(xes[xes <= xdlc[-1]] , yes[xes <= xdlc[-1]] , label = 'Edelstahl', color = 'g')
ax3.plot(xdlc, ydlc, label = 'DLC', color = 'r', linewidth = 0.8)


fig.tight_layout()
fig.savefig('eindringtiefe.pdf', bbox_inches = 'tight', pad_inches = 0)
