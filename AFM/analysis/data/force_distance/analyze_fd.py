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




d_es *= 0.155
d_es_f, d_es_b = np.split(d_es, 2)
F_es_f, F_es_b = np.split(F_es, 2)

fig = plt.figure()

ax1 = fig.add_subplot(131)

ax1.plot(d_es_f, F_es_f, label = 'Vor')
ax1.plot(d_es_b, F_es_b, label = 'Zurück')
ax1.text(0.01 , 0.96, r'(a) Edelstahl', horizontalalignment='left', verticalalignment='top', transform = ax1.transAxes)
ax1.set_ylabel(r'Spannung $ / \si{\volt}$')






d_dlc, F_dlc = np.genfromtxt('force_distance_dlc.csv', unpack = True, delimiter = ';')
d_dlc *= 0.155
d_dlc_f, d_dlc_b = np.split(d_dlc, 2)
F_dlc_f, F_dlc_b = np.split(F_dlc, 2)


ax2 = fig.add_subplot(132, sharey = ax1)
ax2.plot(d_dlc_f, F_dlc_f)
ax2.plot(d_dlc_b, F_dlc_b)
ax2.text(0.01 , 0.96, r'(b) DLC', horizontalalignment='left', verticalalignment='top', transform = ax2.transAxes)
ax2.set_xlabel('Verschiebung $x / \si{\micro\meter}$')







d_t, F_t = np.genfromtxt('force_distance_teflon.csv', unpack = True, delimiter = ';')
d_t *= 0.155
d_t_f, d_t_b = np.split(d_t, 2)
F_t_f, F_t_b = np.split(F_t, 2)



ax3 = fig.add_subplot(133, sharey = ax1)
ax3.plot(d_t_f, F_t_f, label = 'Vor')
ax3.plot(d_t_b, F_t_b, label = 'Zurück')
ax3.text(0.01 , 0.96, r'(c) Teflon', horizontalalignment='left', verticalalignment='top', transform = ax3.transAxes)

ax3.legend(loc = (0.2, 0.6))
fig.tight_layout()
plt.savefig('force_distance.pdf', bbox_inches = 'tight', pad_inches = 0)
