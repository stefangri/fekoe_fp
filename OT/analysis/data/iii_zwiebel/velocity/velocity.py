import numpy as np

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

vehicle_r = np.genfromtxt('vehicel_radius.txt', unpack = True)
print(np.mean(vehicle_r))
quarz_r = np.genfromtxt('quarz_radius.txt', unpack = True)
print(np.mean(quarz_r))


vehicle_r_mum = np.mean(vehicle_r) / quarz_r * 2
print(vehicle_r_mum)


x, y, sum = np.genfromtxt('size_measurement7.dat', unpack = True)
xcal, ycal = np.genfromtxt(f'zcalibration.txt', unpack = True)
x /= xcal
y /= (-ycal)
fs = 10e3
dt = 1 / fs
t = np.arange(0, len(y)*dt, step = dt)
t_in = (1.7, 2.6)
delta_t = t_in[1] - t_in[0]
delta_y = vehicle_r_mum
v = delta_y / delta_t
print(v)


fig = plt.figure(figsize = (5, 2))
ax = fig.add_subplot(111)
ax.plot(t, y)
ax.set_xlabel('$t / \si{\second}$')
ax.set_ylabel('$y / \si{\micro\meter}$')

ax.axvline(t_in[0], color = 'k')
ax.axvline(t_in[1], color = 'k')
fig.savefig('velocity.pdf', bbox_inches = 'tight', pad_inches = 0)
