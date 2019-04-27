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
rcParams['figure.figsize'] =  5.906, 5.906
rcParams['savefig.dpi'] = 300

x1, y1, x2, y2 = np.genfromtxt('kreise_profile.txt', unpack = True)
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

image = mpimg.imread('kreise_250_250_50pps_profile.png')
ax3 = fig.add_subplot(221)
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

image = mpimg.imread('kreise_250_250_100pps_ohne_straingauge.png')
ax4 = fig.add_subplot(222)
ax4.imshow(image)
ax4.set_yticks([])
ax4.set_xticks([])
ax4.text(0.01 , 0.96, r'(b)', horizontalalignment='left', verticalalignment='top', transform = ax4.transAxes, color = 'white')
ax4.set_xticks([0, 250, 500])
ax4.set_yticks([0, 250, 500])
ax4.set_xticklabels([0, 10, 20])
ax4.set_yticklabels([20, 10, 0])
ax4.set_xlabel('$x / \si{\micro\meter}$')
ax4.set_ylabel('$y / \si{\micro\meter}$')

ax1 = fig.add_subplot(223)

ax1.plot(x1, y1, label = 'Profil horizontal')
ax1.set_xlim(x1[0], x1[-8])
ax1.set_xlabel(r'Position $x / \si{\micro\meter}$ ')
ax1.set_ylabel('Höhe $y / \si{\micro\meter}$')
ax1.legend(loc = 'lower right')
ax1.text(0.01 , 0.96, r'(b)', horizontalalignment='left', verticalalignment='top', transform = ax1.transAxes)

featurex, featurey = np.genfromtxt('features1_kreise_profile.txt', unpack = True)
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

yplot = np.linspace(3.95, 4.3, 100)
#ax1.set_ylim(3.95, 4.3)
for feature in featurex:
    ax1.axvline(x = feature, linestyle = '--', color = 'black', linewidth = 0.5)
    ax1.fill_betweenx(yplot, feature - 0.2, feature + 0.2, color = 'grey', alpha = 0.5, linewidth = 0)

upper_mask = (x1 > featurex[-1] + 0.2)
lower_mask = x1 < featurex[0] - 0.2
for i in [0, 2, 4]:
    upper_mask = upper_mask | ((x1 > featurex[i] + 0.2) & (x1 < featurex[i + 1] - 0.2))

for i in [1, 3, 5]:
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
print(np.mean(lin(x1, *upper_params) - lin(x1, *lower_params)))


ax1.legend(loc = 'lower right')



















ax2 = fig.add_subplot(224, sharey = ax1)

ax2.plot(x2, y2, label = 'Profil vertikal')
ax2.set_xlim(x2[0], x2[-8])
ax2.set_xlabel(r'Position $x / \si{\micro\meter}$ ')
ax2.legend(loc = 'lower right')
ax2.text(0.01 , 0.96, r'(c)', horizontalalignment='left', verticalalignment='top', transform = ax2.transAxes)
#ax2.set_ylabel('Höhe $y / \si{\micro\meter}$')



featurex, featurey = np.genfromtxt('features2_kreise_profile.txt', unpack = True)
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

yplot = np.linspace(3.95, 4.3, 100)
ax2.set_ylim(3.95, 4.3)
for feature in featurex:
    ax2.axvline(x = feature, linestyle = '--', color = 'black', linewidth = 0.5)
    ax2.fill_betweenx(yplot, feature - 0.4, feature + 0.4, color = 'grey', alpha = 0.5, linewidth = 0)

upper_mask =  (x2 > featurex[-1] + 0.4)
lower_mask = x2 == y2
for i in [1, 3, 5]:
    upper_mask = upper_mask | ((x2 > featurex[i] + 0.4) & (x2 < featurex[i + 1] - 0.4))

for i in [0, 2, 4, 6]:
    lower_mask = lower_mask | ((x2 > featurex[i] + 0.4) & (x2 < featurex[i + 1] - 0.4))

def lin(x, a, b):
    return a * x + b

upper_params, upper_cov = curve_fit(lin, x2[upper_mask], y2[upper_mask])
lower_params, lower_cov = curve_fit(lin, x2[lower_mask], y2[lower_mask])

ax2.plot(x2, lin(x2, *upper_params), color = 'black', linestyle = '-', linewidth = 1, label = 'Regression')
ax2.plot(x2, lin(x2, *lower_params), color = 'black', linestyle = '-', linewidth = 1)

upper_params = correlated_values(upper_params, upper_cov)
lower_params = correlated_values(lower_params, lower_cov)

print(upper_params)
print(lower_params)
print(np.mean(lin(x2, *upper_params) - lin(x2, *lower_params)))


ax2.legend(loc = 'lower left')


fig.tight_layout()
fig.savefig('kreise_profile.pdf', bbox_inches = 'tight', pad_inches = 0)
