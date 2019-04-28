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
rcParams['figure.figsize'] =  5.906, 2.5
rcParams['savefig.dpi'] = 300

x1, y1 = np.genfromtxt('cd_profile.txt', unpack = True)
x1 *= 1e6
y1 *= 1e6

fig = plt.figure()

image = mpimg.imread('cd_10_10_100pps_profiles.png')
ax3 = fig.add_subplot(131)
ax3.imshow(image)
ax3.set_yticks([])
ax3.set_xticks([])
ax3.text(0.01 , 0.96, r'(a) CD', horizontalalignment='left', verticalalignment='top', transform = ax3.transAxes, color = 'black')
ax3.set_xticks([0, 2000, 4000])
ax3.set_yticks([0, 2000, 4000])
ax3.set_xticklabels([0, 5, 10])
ax3.set_yticklabels([10, 5, 0])
ax3.set_xlabel('$x / \si{\micro\meter}$')
ax3.set_ylabel('$y / \si{\micro\meter}$')




image = mpimg.imread('../dvd/dvd_250_250_100pps.png')
ax2 = fig.add_subplot(133)
ax2.imshow(image)
ax2.set_yticks([])
ax2.set_xticks([])
ax2.text(0.01 , 0.96, r'(c) DVD', horizontalalignment='left', verticalalignment='top', transform = ax2.transAxes, color = 'black')
ax2.set_xticks([0, 250, 500])
ax2.set_yticks([0, 250, 500])
ax2.set_xticklabels([0, 2.5, 5])
ax2.set_yticklabels([5, 2.5, 0])
ax2.set_xlabel('$x / \si{\micro\meter}$')
ax2.set_ylabel('$y / \si{\micro\meter}$')

ax1 = fig.add_subplot(132)

ax1.plot(x1, y1, label = 'Profil')
ax1.set_xlim(x1[0], x1[-1])
ax1.set_xlabel(r'Position $x / \si{\micro\meter}$ ')
ax1.set_ylabel('Höhe $y / \si{\micro\meter}$')

ax1.text(0.05 , 0.96, r'(b)', horizontalalignment='left', verticalalignment='top', transform = ax1.transAxes)


featurex, featurey = np.genfromtxt('features_cd_profile.txt', unpack = True)
featurex *= 1e6



yplot = np.linspace(ax1.get_ylim()[0], ax1.get_ylim()[-1], 100)
ax1.set_ylim(yplot[0], yplot[-1])
for feature in featurex:
    ax1.axvline(x = feature, linestyle = '--', color = 'black', linewidth = 0.5)
    ax1.fill_betweenx(yplot, feature - 0.2, feature + 0.2, color = 'grey', alpha = 0.5, linewidth = 0)

upper_mask = (x1 < featurex[0] - 0.2) | (x1 > featurex[-1] + 0.2)
lower_mask = x1 == y1
for i in [1, 3]:
    upper_mask = upper_mask | ((x1 > featurex[i] + 0.2) & (x1 < featurex[i + 1] - 0.2))

for i in [0, 2, 4]:
    lower_mask = lower_mask | ((x1 > featurex[i] + 0.2) & (x1 < featurex[i + 1] - 0.2))

def lin(x, a, b):
    return a * x + b

upper_params, upper_cov = curve_fit(lin, x1[upper_mask], y1[upper_mask])
lower_params, lower_cov = curve_fit(lin, x1[lower_mask], y1[lower_mask])



ax1.plot(x1, lin(x1, *upper_params), color = 'black', linestyle = '-', linewidth = 1)
ax1.plot(x1, lin(x1, *lower_params), color = 'black', linestyle = '-', linewidth = 1)

upper_params = correlated_values(upper_params, upper_cov)
lower_params = correlated_values(lower_params, lower_cov)

print(upper_params)
print(lower_params)
print(upper_params[1] - lower_params[1])



ax1.legend(loc = 'upper right')
fig.tight_layout()
fig.savefig('cd_profil.pdf', bbox_inches = 'tight', pad_inches = 0)
