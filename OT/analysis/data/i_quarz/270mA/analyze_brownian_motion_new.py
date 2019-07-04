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



def power_density(f, A, f0):
    return A / (f**2 / f0 + f0)

name = os.getcwd().split('/')[-1]
with open(f'results/results_without_force_{name}.txt', 'w') as file:
    file.write(f'#Amplitudes and roll-off frequencies for laser power {name} \n')

#ANALYZE DATA WITHOUT FORCE
x, y, sum = np.genfromtxt(f'{name}_ohne.dat', unpack = True)
xcal, ycal = np.genfromtxt(f'zcalibration.txt', unpack = True)
x /= xcal
y /= ycal
fs = 10e3

f, Pxx_den = signal.periodogram(x - np.mean(x), fs)
dt = 1 / fs
t = np.arange(0, len(x)*dt, step = dt)

fig = plt.figure(figsize = (5.7, 7.8))
ax1 = fig.add_subplot(411)
ax1.plot(t, x, label = 'Daten')
ax1.set_xlabel(r'$t / \si{\second}$')
ax1.set_ylabel(r'$x / \si{\micro\meter}$')
ax1.set_xlim(t[0], t[-1])
ax1.legend(loc = 'lower left')
ax1.text(0.03 , 0.96, r'(a)', horizontalalignment='left', verticalalignment='top', transform = ax1.transAxes, color = 'black')
#fig.savefig(f'xplot_{name}.pdf', bbox_inches = 'tight', pad_inches = 0)


ax2 = fig.add_subplot(412)

ffcal, other, xfcal, yfcal = np.genfromtxt(f'force_calibration_{name}.FDdat', unpack = True)


fit_lim = (10, 3000)
start_params = [1e-4, 2]
ax2.plot(f, Pxx_den)
Pxx_den[f < 3] = Pxx_den[f>3][0]
params, cov = curve_fit(power_density, f[(f > fit_lim[0]) & (f < fit_lim[1])], Pxx_den[(f > fit_lim[0]) & (f < fit_lim[1])], p0 = start_params)
Ax = ufloat(params[0], np.sqrt(np.diag(cov))[0])
f0x = Q_(ufloat(params[1], np.sqrt(np.diag(cov))[1]), 'hertz')
with open(f'results/results_with_force_{name}.txt', 'a') as file:
    file.write(f'Ax \t {params[0]} \t {np.sqrt(np.diag(cov))[0]} \n')
    file.write(f'f0x \t {params[1]} \t {np.sqrt(np.diag(cov))[1]} \n')

fplot = np.logspace(-1, 5)

#params2, cov2 = curve_fit(power_density, x[(x > 10) & (x < 800)], z[(x > 10) & (x < 800)], p0 = [1e-8, 10])
#print(params2)


#ax2.plot(ffcal, xfcal, label = 'Software Daten')
#ax.plot(fplot, power_density(fplot, *params2), 'b-')
ax2.plot(fplot, power_density(fplot, *params), 'r-', label = 'Fit')
ax2.set_xlabel(r'$f / \si{\hertz}$')
ax2.set_ylabel(r'PSD$_x / \si{\micro\meter\squared \second}$')
ax2.set_xscale('log')
ax2.set_yscale('log')
ylim = [1e-15, 1e-2]
ax2.set_xlim(f[0], f[-1])
ax2.set_ylim(ylim)
ax2.legend(loc = 'lower left')
ax2.fill_betweenx(y = np.linspace(*ylim), x1 = fit_lim[0], x2 = fit_lim[1], color = 'grey', alpha = 0.5)
ax2.text(0.01 , 0.96, r'(b)', horizontalalignment='left', verticalalignment='top', transform = ax2.transAxes, color = 'black')
#ax2.savefig(f'xpowerdensity_{name}.png')



f, Pxx_den = signal.periodogram(y - np.mean(y), fs)
dt = 1 / fs
t = np.arange(0, len(y)*dt, step = dt)

ax3 = fig.add_subplot(413)

ax3.plot(t, y, color = colors[2], label = 'Daten')
ax3.set_xlim(t[0], t[-1])
ax3.set_xlabel(r'$t / \si{\second}$')
ax3.set_ylabel(r'$y / \si{\micro\meter}$')
ax3.legend(loc = 'lower left')
ax3.text(0.01 , 0.96, r'(c)', horizontalalignment='left', verticalalignment='top', transform = ax3.transAxes, color = 'black')




ax4 = fig.add_subplot(414)
ffcal, other, xfcal, yfcal = np.genfromtxt(f'force_calibration_{name}.FDdat', unpack = True)


fit_lim = (100, 3000)
start_params = [1e-4, 2]
ax4.plot(f, Pxx_den, color = colors[2])
Pxx_den[f < 3] = Pxx_den[f>3][0]
params, cov = curve_fit(power_density, f[(f > fit_lim[0]) & (f < fit_lim[1])], Pxx_den[(f > fit_lim[0]) & (f < fit_lim[1])], p0 = start_params)
Ay = ufloat(params[0], np.sqrt(np.diag(cov))[0])
f0y = Q_(ufloat(params[1], np.sqrt(np.diag(cov))[1]), 'hertz')
with open(f'results/results_with_force_{name}.txt', 'a') as file:
    file.write(f'Ay \t {params[0]} \t {np.sqrt(np.diag(cov))[0]} \n')
    file.write(f'f0y \t {abs(params[1])} \t {np.sqrt(np.diag(cov))[1]} \n')


fplot = np.logspace(-1, 5)

#ax4.plot(ffcal, xfcal)
#ax.plot(fplot, power_density(fplot, *params2), 'b-')
ax4.plot(fplot, power_density(fplot, *params), 'r-', label = 'Fit')
ax4.set_xlabel(r'$f / \si{\hertz}$')
ax4.set_ylabel(r'PSD$_y/ \si{\micro\meter\squared \second}$')
ax4.set_xscale('log')
ax4.set_yscale('log')
ylim = [1e-15, 1e-2]
ax4.set_ylim(ylim)
ax4.set_xlim(f[0], f[-1])
ax4.fill_betweenx(y = np.linspace(*ylim), x1 = fit_lim[0], x2 = fit_lim[1], color = 'grey', alpha = 0.5)
ax4.legend(loc = 'lower left')
ax4.text(0.01 , 0.9, r'(d)', horizontalalignment='left', verticalalignment='top', transform = ax4.transAxes, color = 'black')
#ax.savefig(f'ypowerdensity_{name}.png')

fig.tight_layout()






# Calculate k_b and spring constant



d = Q_(2, 'micrometer')
eta = Q_(8.9e-4, 'pascal * second')
T = Q_(20, 'celsius').to('kelvin')


beta = 3 * np.pi * eta * d
kx = (2 * np.pi * beta * f0x).to('newton / meter')

#kb = Q_(1.4e-23, 'joule / kelvin')#3.204973953607321e-26
#kx_alternative = (2 * kb * T / np.pi / Q_(Ax, 'micrometer**2')).to('newton / meter')
#print(kx_alternative)

ky = abs((2 * np.pi * beta * f0y)).to('newton / meter')


#load data from measurement without force
x, y, sum = np.genfromtxt(f'{name}_ohne.dat', unpack = True)
x /= xcal
y /= ycal
dt = 1 / fs
t = np.arange(0, len(x)*dt, step = dt)



fig2 = plt.figure(figsize = (5.7, 2))
#ax5 = fig2.add_subplot(221)
#ax5.plot(t, x - np.mean(x))
#ax5.set_xlim(t[0], t[-1])
#ax5.set_xlabel(r'$t / \si{\second}$')
#ax5.set_ylabel(r'$x - \langle x \rangle / \si{\micro\meter}$')
#ax5.text(0.01 , 0.96, r'(a)', horizontalalignment='left', verticalalignment='top', transform = ax5.transAxes, color = 'black')
#
#ax6 = fig2.add_subplot(223)
#ax6.plot(t, y - np.mean(y), color = colors[1])
#ax6.set_xlim(t[0], t[-1])
#ax6.set_xlabel(r'$t / \si{\second}$')
#ax6.set_ylabel(r'$y - \langle y \rangle / \si{\micro\meter}$')
#ax6.text(0.01 , 0.96, r'(c)', horizontalalignment='left', verticalalignment='top', transform = ax6.transAxes, color = 'black')


ax7 = fig2.add_subplot(121)
ax7.hist(x - np.mean(x), bins = 30, histtype = 'step')
ax7.set_xlabel(r'$x - \langle x \rangle / \si{\micro\meter}$')
ax7.text(0.01 , 0.96, r'(a)', horizontalalignment='left', verticalalignment='top', transform = ax7.transAxes, color = 'black')
ax7.set_ylabel('Häufigkeit')

ax8 = fig2.add_subplot(122)
ax8.hist(y - np.mean(y), bins = 30, histtype = 'step', color = colors[1])
ax8.set_xlabel(r'$y - \langle y \rangle / \si{\micro\meter}$')
ax8.text(0.01 , 0.96, r'(b)', horizontalalignment='left', verticalalignment='top', transform = ax8.transAxes, color = 'black')
ax8.set_ylabel('Häufigkeit')
fig2.tight_layout()






varx = Q_(np.var(x - np.mean(x)), 'micrometer**2')
vary = Q_(np.var(y - np.mean(y)), 'micrometer**2')
print(varx, vary)



kb_x = (kx * varx / T).to('joule / kelvin')
kb_y = (ky * vary / T).to('joule / kelvin')


with open(f'results/results_spring_constant_{name}.txt', 'w') as file:
    file.write(f'#Spring constants for laser power {name} \n')

with open(f'results/results_spring_constant_{name}.txt', 'a') as file:
    file.write(f'kx/N/m \t {kx.magnitude.n} \t {kx.magnitude.s} \n')
    file.write(f'ky/N/m \t {ky.magnitude.n} \t {ky.magnitude.s} \n')


with open(f'results/results_boltzmann_constant_{name}.txt', 'w') as file:
    file.write(f'#boltzmann constant for laser power {name} \n')

with open(f'results/results_boltzmann_constant_{name}.txt', 'a') as file:
    file.write(f'kbx/J/K \t {kb_x.magnitude.n} \t {kb_x.magnitude.s} \n')
    file.write(f'kby/J/K \t {kb_y.magnitude.n} \t {kb_y.magnitude.s} \n')



# DATA WITH FORCE
x, y, sum = np.genfromtxt(f'{name}_x.dat', unpack = True)
x /= xcal
y /= ycal
t = np.arange(0, len(x)*dt, step = dt)

v = Q_(4, 'micrometer / second')
kx_force = np.mean(beta * v / Q_(x- np.mean(x), 'micrometer')).to('newton / meter')

fig3 = plt.figure(figsize = (5.7, 3.9))
ax9 = fig3.add_subplot(211)
ax9.plot(t, x, label = 'Daten')
ax9.set_xlabel(r'$t / \si{\second}$')
ax9.set_ylabel(r'$x / \si{\micro\meter}$')
ax9.set_xlim(t[0], t[-1])
ax9.legend(loc = 'lower left')
ax9.text(0.03 , 0.96, r'(a)', horizontalalignment='left', verticalalignment='top', transform = ax9.transAxes, color = 'black')



x, y, sum = np.genfromtxt(f'{name}_y.dat', unpack = True)
x /= xcal
y /= ycal
t = np.arange(0, len(y)*dt, step = dt)


ky_force = np.mean(beta * v / Q_(y- np.mean(y), 'micrometer')).to('newton / meter')


ax10 = fig3.add_subplot(212)
ax10.plot(t, y, label = 'Daten', color = colors[1])
ax10.set_xlabel(r'$t / \si{\second}$')
ax10.set_ylabel(r'$y / \si{\micro\meter}$')
ax10.set_xlim(t[0], t[-1])
ax10.legend(loc = 'lower left')
ax10.text(0.03 , 0.96, r'(b)', horizontalalignment='left', verticalalignment='top', transform = ax10.transAxes, color = 'black')

kb_xforce = (kx_force * varx / T).to('joule / kelvin')
kb_yforce = (ky_force * vary / T).to('joule / kelvin')

with open(f'results/results_boltzmann_constant_{name}.txt', 'a') as file:
    file.write(f'kbx/J/Kforce \t {kb_xforce.magnitude} \t {0} \n')
    file.write(f'kby/J/Kforce \t {kb_yforce.magnitude} \t {0} \n')   


with open(f'results/results_spring_constant_{name}.txt', 'a') as file:
    file.write(f'kx/N/mwithforce \t {kx_force.magnitude} \t {0} \n')
    file.write(f'ky/N/mwithforce \t {ky_force.magnitude} \t {0} \n')
















fig.savefig(f'results/without_force_{name}.pdf', bbox_inches = 'tight', pad_inches = 0)
fig2.savefig(f'results/without_force_histogram_{name}.pdf', bbox_inches = 'tight', pad_inches = 0)
fig3.savefig(f'results/with_force_{name}.pdf', bbox_inches = 'tight', pad_inches = 0)
