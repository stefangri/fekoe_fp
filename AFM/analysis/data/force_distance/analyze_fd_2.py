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
rcParams['figure.figsize'] =  5.906, 3
rcParams['savefig.dpi'] = 300


d_es, F_es = np.genfromtxt('force_distance_edelstahl.csv', unpack = True, delimiter = ';')
d_es *= (20 / 75)
d_es_f, d_es_b = np.split(d_es, 2)
F_es_f, F_es_b = np.split(F_es, 2)

features = np.genfromtxt('features_force_distance_edelstahl.txt')
pulloff = features[1]
snapin = features[0]
pulloff[0] *= (20 / 75)
snapin[0] *= (20 / 75)
zero, cov = curve_fit(lambda x, c: [c for i in x], d_es[d_es < 0.5], F_es[d_es < 0.5])
lin = lambda x, a, b: a * x + b
lin_forward, lin_cov = curve_fit(lin, d_es_f[d_es_f > snapin[0]], F_es_f[d_es_f > snapin[0]])
lin_backward, lin_cov = curve_fit(lin, d_es_b[d_es_b > pulloff[0]], F_es_b[d_es_b > pulloff[0]])

np.savetxt('edelstahl_repulsive.txt', np.column_stack([d_es_f[d_es_f > snapin[0]], F_es_f[d_es_f > snapin[0]]]))


fig = plt.figure()

ax1 = fig.add_subplot(311)

xplot = np.linspace(0.5, 3, 100)
ax1.plot(*pulloff, 'ro', label = 'Pull off')
ax1.plot(*snapin, 'go', label = 'Snap in')
ax1.axhline(y = zero, color = 'grey', linestyle = '--', label = 'Ruhelage')
ax1.plot(xplot, lin(xplot, *lin_forward), color = 'black', label = 'Regression')
ax1.plot(xplot, lin(xplot, *lin_backward), color = 'black', linestyle = '-')

ax1.set_ylim(-0.5, 0.8)
ax1.plot(d_es_f, F_es_f, label = 'Vor')
ax1.plot(d_es_b, F_es_b, label = 'Zurück')
ax1.text(0.01 , 0.96, r'(a) Edelstahl', horizontalalignment='left', verticalalignment='top', transform = ax1.transAxes)
#ax1.set_ylabel(r'Kraft $F / \si{\newton}$')
ax1.legend(loc = 'lower right')






d_dlc, F_dlc = np.genfromtxt('force_distance_dlc.csv', unpack = True, delimiter = ';')
d_dlc *= (20 / 75)
d_dlc_f, d_dlc_b = np.split(d_dlc, 2)
F_dlc_f, F_dlc_b = np.split(F_dlc, 2)

features = np.genfromtxt('features_force_distance_dlc.txt')
pulloff = features[1]
snapin = features[0]
pulloff[0] *= (20 / 75)
snapin[0] *= (20 / 75)
np.savetxt('dlc_repulsive.txt', np.column_stack([d_dlc_f[d_dlc_f > snapin[0]], F_dlc_f[d_dlc_f > snapin[0]]]))

ax2 = fig.add_subplot(312, sharex = ax1)
ax2.plot(d_dlc_f, F_dlc_f)
ax2.plot(d_dlc_b, F_dlc_b)
ax2.text(0.01 , 0.96, r'(b) DLC', horizontalalignment='left', verticalalignment='top', transform = ax2.transAxes)
ax2.set_ylabel(r'Spannung $/ \si{\volt}$')






d_t, F_t = np.genfromtxt('force_distance_teflon.csv', unpack = True, delimiter = ';')
d_t *= (20 / 75)
d_t_f, d_t_b = np.split(d_t, 2)
F_t_f, F_t_b = np.split(F_t, 2)


features = np.genfromtxt('features_force_distance_teflon.txt')
pulloff = features[1]
snapin = features[0]
pulloff[0] *= (20 / 75)
snapin[0] *= (20 / 75)
zero, cov = curve_fit(lambda x, c: [c for i in x], d_t[d_t < 0.5], F_t[d_t < 0.5])
np.savetxt('teflon_repulsive.txt', np.column_stack([d_t_f[d_t_f > snapin[0]], F_t_f[d_t_f > snapin[0]]]))

fig = plt.figure()
ax3 = fig.add_subplot(111)

ax3.axhline(y = zero, color = 'grey', linestyle = '--', label = 'Ruhelage')
ax3.plot(*pulloff, 'ro', label = 'Pull off')
ax3.plot(*snapin, 'go', label = 'Snap in')

ax3.plot(d_t_f, F_t_f, label = 'Vor')
ax3.plot(d_t_b, F_t_b, label = 'Zurück')
#ax3.text(0.01 , 0.96, r'(c) Teflon', horizontalalignment='left', verticalalignment='top', transform = ax3.transAxes)
ax3.set_xlabel('Verschiebung $x / \si{\micro\meter}$')
ax3.set_ylabel(r'Spannung $/ \si{\volt}$')
ax3.text(3.2, -0.2, 'Kontakt')

ax3.legend()
fig.tight_layout()
plt.savefig('force_distance_teflon.pdf', bbox_inches = 'tight', pad_inches = 0)
