import numpy as np
from tab2tex import *
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




#get calibration data
slopex = []
slopey = []
for power in powers:
    x, y = np.genfromtxt(f'{power}mA/zcalibration.txt', unpack = True)
    slopex.append(x)
    slopey.append(y)
slopex = np.array(slopex)
slopey = np.array(slopey)
make_table(header = [r'$I$ / \milli\ampere', r'$s_x$ / \volt\per\micro\meter', r'$s_y$ / \volt\per\micro\meter'],
places = [3.0, 1.3, 1.3], data = [powers, slopex, slopey], caption = 'Umrechnungsfaktoren zwischen Spannung der Viersegmentdiode und Strecke für die verwendeten Pumpströme.',
label = 'tab: zcalibration', filename = 'zcalibration_tab.tex')

#mean of boltzmann constant
kbx = []
kby = []
kbxforce = []
kbyforce = []
for power in powers:
    foo, kb, baz = np.genfromtxt(f'{power}mA/results/results_boltzmann_constant_{power}mA.txt', unpack = True)
    kbx.append(kb[0])
    kby.append(kb[1])
    kbxforce.append(kb[2])
    kbyforce.append(kb[3])
kbx = np.array(kbx)
kby = np.array(kby)
kbxforce = np.array(kbxforce)
kbyforce = np.array(kbyforce)

print('kbmean', np.mean([kbx, kby]))
print('kbmeanforce', np.mean([kbxforce, kbyforce]))



####################################

kx = []
ky = []
kxforce = []
kyforce = []

for power in powers:
    na, k, kerr = np.genfromtxt(f'{power}mA/results/results_spring_constant_{power}mA.txt', unpack = True)
    kx.append(k[0])
    ky.append(k[1])
    kxforce.append(k[2])
    kyforce.append(k[3])
kx = np.array(kx)
ky = np.array(ky)
kxforce = np.array(kxforce)
kyforce = np.array(kyforce)

paramsx, covx = curve_fit(lin, powers, kx * 1e7)
paramsy, covy = curve_fit(lin, powers, ky * 1e7)

print(paramsx, paramsy)

make_table(header = [r'$I$ / \milli\ampere', r'$k_x$ / 10^{-7}\newton\per\meter', r'$k_y$ / 10^{-7}\newton\per\meter', r'$k_{x, F}$ / 10^{-7}\newton\per\meter', r'$k_{y, F}$ / 10^{-7}\newton\per\meter'],
places = [3.0, 3.1, 3.1, 3.1, 3.1], data = [powers, kx*1e7, ky*1e7, abs(kxforce)*1e7, abs(kyforce)*1e7], caption = 'Ergebnisse für die Fallensteifigkeit in $x$- und $y$-Richtung. Die Ergebnisse wurden aus den Messungen ohne Krafteinwirkung sowie aus den Messungen mit Krafteinwirkung (Index F) gewonnen.',
label = 'tab: quarz_result', filename = 'quarz_result_tab.tex')



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



fig = plt.figure(figsize = (5, 2))
ax = fig.add_subplot(111)
ax.plot(powers, abs(kxforce) * 1e7, 'o', fillstyle = 'none', label = 'Daten $k_x$')
ax.plot(powers, abs(kyforce) * 1e7, 'o', fillstyle = 'none', label = 'Daten $k_y$')
ax.set_ylabel(r'$k / \si{\num{e-7} \newton \per \meter}$')
ax.set_xlabel(r'$I / \si{\milli\ampere}$')
ax.set_xticks(powers)
ax.set_xlim(50, 490)
ax.legend()
fig.savefig('k_power_series_force.pdf', bbox_inches = 'tight', pad_inches = 0)
