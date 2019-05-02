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

x1, y1 = np.genfromtxt('linien_profil.txt', unpack = True)
x1 *= 1e6
y1 *= 1e6

fig = plt.figure()

image = mpimg.imread('linien_250_250_100pps_profile.png')
ax3 = fig.add_subplot(121)
ax3.imshow(image)
ax3.set_yticks([])
ax3.set_xticks([])
ax3.text(0.01 , 0.96, r'(a)', horizontalalignment='left', verticalalignment='top', transform = ax3.transAxes, color = 'white')
ax3.set_xticks([0, 2000, 4000])
ax3.set_yticks([0, 2000, 4000])
ax3.set_xticklabels([0, 10, 20])
ax3.set_yticklabels([20, 10, 0])
ax3.set_xlabel('$x / \si{\micro\meter}$')
ax3.set_ylabel('$y / \si{\micro\meter}$')

ax1 = fig.add_subplot(122)




ax1.plot(x1, y1, label = 'Profil')
ax1.set_xlim(x1[0], x1[-1])
ax1.set_xlabel(r'Position $x / \si{\micro\meter}$ ')
ax1.set_ylabel('Höhe $y / \si{\micro\meter}$')

ax1.text(0.05 , 0.96, r'(b)', horizontalalignment='left', verticalalignment='top', transform = ax1.transAxes)


featurex, featurey = np.genfromtxt('features_linien_profil.txt', unpack = True)
featurex *= 1e6

horizontal_dist = np.diff(featurex)
gaps = horizontal_dist[::2]
steps = horizontal_dist[1::2]
mean_gap = ufloat(np.mean(gaps), np.std(gaps))
mean_step = ufloat(np.mean(steps), np.std(steps))
structure_distance = mean_gap + mean_step
print(mean_gap)
print(mean_step)
print(structure_distance)


yplot = np.linspace(ax1.get_ylim()[0], ax1.get_ylim()[-1], 100)
ax1.set_ylim(yplot[0], yplot[-1])
for feature in featurex:
    ax1.axvline(x = feature, linestyle = '--', color = 'black', linewidth = 0.5)
    ax1.fill_betweenx(yplot, feature - 0.2, feature + 0.2, color = 'grey', alpha = 0.5, linewidth = 0)

upper_mask = (x1 < featurex[0] - 0.2) | (x1 > featurex[-1] + 0.2)
lower_mask = x1 == y1
for i in [1, 3, 5]:
    upper_mask = upper_mask | ((x1 > featurex[i] + 0.2) & (x1 < featurex[i + 1] - 0.2))

for i in [0, 2, 4, 6]:
    lower_mask = lower_mask | ((x1 > featurex[i] + 0.2) & (x1 < featurex[i + 1] - 0.2))

def lin(x, a, b):
    return a * x + b

upper_params, upper_cov = curve_fit(lin, x1[upper_mask], y1[upper_mask])
lower_params, lower_cov = curve_fit(lin, x1[lower_mask], y1[lower_mask])



ax1.plot(x1, lin(x1, *upper_params), color = 'black', linestyle = '-', linewidth = 1, label = 'Regression')
ax1.plot(x1, lin(x1, *lower_params), color = 'black', linestyle = '-', linewidth = 1)

upper_params = correlated_values(upper_params, upper_cov)
lower_params = correlated_values(lower_params, lower_cov)

print(upper_params)
print(lower_params)
print(upper_params[1] - lower_params[1])



ax1.legend(loc = 'upper right')
fig.tight_layout()
fig.savefig('linien_profil.pdf', bbox_inches = 'tight', pad_inches = 0)
