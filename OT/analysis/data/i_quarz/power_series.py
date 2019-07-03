import numpy as np

from scipy import signal
from scipy.optimize import curve_fit
import os
from uncertainties import ufloat
from pint import UnitRegistry
Q_ = UnitRegistry().Quantity

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
rcParams['axes.formatter.use_mathtext'] = True
rcParams['legend.fontsize'] = 10
rcParams['savefig.dpi'] = 300
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

def lin(x, a, b):
    return a * x + b

powers = np.array([70, 170, 270, 370, 470])


kx = []
ky = []

for power in powers:
    na, k, kerr = np.genfromtxt(f'{power}mA/results/results_spring_constant_{power}mA.txt', unpack = True)
    kx.append(k[0])
    ky.append(k[1])
kx = np.array(kx)
ky = np.array(ky)

paramsx, covx = curve_fit(lin, powers, kx* 1e7)
paramsy, covy = curve_fit(lin, powers, ky* 1e7)







pplot = np.linspace(50, 490, 100)

fig = plt.figure(figsize = (5, 3))
ax = fig.add_subplot(111)
ax.plot(powers, kx * 1e7, 'o', fillstyle = 'none', label = 'Daten $k_x$')
ax.plot(pplot, lin(pplot, *paramsx), color = colors[0], linestyle = '--', label = 'Fit $k_x$')
ax.plot(powers, ky * 1e7, 'o', fillstyle = 'none', label = 'Daten $k_y$')
ax.plot(pplot, lin(pplot, *paramsy), color = colors[1], linestyle = '--', label = 'Fit $k_y$')
ax.set_ylabel(r'$k / \si{\num{e-7} \newton \per \meter}$')
ax.set_xlabel(r'$I / \si{\milli\ampere}$')
ax.set_xticks(powers)
ax.set_xlim(50, 490)
ax.legend()
fig.savefig('k_power_series.pdf', bbox_inches = 'tight', pad_inches = 0)
