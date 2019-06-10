import numpy as np
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as stds
from uncertainties import correlated_values
import math
from scipy.optimize import curve_fit
from pint import UnitRegistry
import latex as l
u = UnitRegistry()
Q_ = u.Quantity
import os

def abs_path(filename):
    return os.path.join(os.path.dirname(__file__), filename)

import matplotlib as mlp
mlp.use('pgf')
import matplotlib.pyplot as plt
from matplotlib import rcParams

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
rcParams['savefig.dpi'] = 300


import numpy as np
import matplotlib.pyplot as plt
#from scipy.special import riccati_yn
#from scipy.special import riccati_jn
from scipy.special import spherical_jn
from scipy.special import spherical_yn
from scipy.optimize import curve_fit
from scipy import interpolate


def riccati_jn(n, x, derivative = False):
    if derivative == True:
        return spherical_jn(n, x, derivative = False) + x * spherical_jn(n, x, derivative = True)
    else:
        return x * spherical_jn(n, x, derivative = False)


def riccati_yn(n, x, derivative = False):
    if derivative == True:
        return (spherical_jn(n, x, derivative = False) + (0 +1j) * spherical_yn(n, x, derivative = False)) + x * (spherical_jn(n, x, derivative = True) + (0 +1j) * spherical_yn(n, x, derivative = True))
    else:
        return x * (spherical_jn(n, x, derivative = False) + (0 +1j) * spherical_yn(n, x, derivative = False))



xlight, ylight = np.genfromtxt(abs_path('data/light_source_spectrum.txt'), unpack = True)
ylight /= max(ylight) #norm spectrum
light_inter = interpolate.interp1d(xlight, ylight, kind = 'linear') #interpolate spectrum




l, ref_n, k = np.genfromtxt(abs_path('data/refrective_index.txt'), unpack = True) #data for refrective index
l *= 1e3 #convert lambda to nanometer

lam_plot = np.linspace(400, 800, 1000) #wavelength array for plots


n_inter = interpolate.interp1d(l, ref_n, kind = 'quadratic') #interpolate real part of refrective index
#plt.plot(lam_plot, n_inter(lam_plot))
#plt.plot(l, ref_n, 'ro')
#plt.show()


k_inter = interpolate.interp1d(l, k, kind = 'quadratic') #interpolate imaginary part of refrective index
#plt.plot(lam_plot, k_inter(lam_plot))
#plt.plot(l, k, 'ro')
#plt.show()

x104, y104 = np.genfromtxt(abs_path('data/sag5b_0001.txt'), unpack = True) #load spectrum data
x14, y14 = np.genfromtxt(abs_path('data/sag6_0001.txt'), unpack = True)
xdark, ydark = np.genfromtxt(abs_path('data/sag7_0001.txt'), unpack = True)



N_medium = 1.5 #refrective index of glas (assume k = 0)



def x(lam, size): #wavelegth dependent parameter for mie scattering function, size = diameter of spherical particle
    return  np.pi / lam  * N_medium * size



def a_n(x, m, n): #first mie coefficient
    return np.array([(m * riccati_jn(n, m * x) * riccati_jn(n, x, True) - riccati_jn(n, x) * riccati_jn(n, m * x, True)) / (m * riccati_jn(n, m * x) * riccati_yn(n, x, True) - riccati_jn(n, m * x, True) * riccati_yn(n, x)) for x, m in zip(x, m)])


def b_n(x, m, n): #second mie coefficient
    return np.array([(riccati_jn(n, m * x) * riccati_jn(n, x, True) - m * riccati_jn(n, x) * riccati_jn(n, m * x, True)) / (riccati_jn(n, m * x) * riccati_yn(n, x, True) - m * riccati_jn(n, m * x, True) * riccati_yn(n, x)) for x, m in zip(x, m)])


def Q_sc(lam, A, size): #scattering cross section
    n = 1 #highest order of mie row
    N_partikel = n_inter(lam) * (1 + 0j) + k_inter(lam)*(0 + 1j)
    m = N_partikel / N_medium
    x_value = x(lam, size)
    a = np.array([a_n(x_value, m, n) for n in range(1, n + 1)])
    b = np.array([b_n(x_value, m, n) for n in range(1, n + 1)])
    return A * 2 / x_value**2 * np.sum( (2 * n + 1) * np.real(a * np.conj(a) + b * np.conj(b)), axis = 0)



def sum_light_mie(lam, A, size, B): #define sum of light source spectrum and mie cross section as model for all data
    return Q_sc(lam, A, size) + B * light_inter(lam)


params1, cov = curve_fit(sum_light_mie, x104[::20], (y104 - ydark)[::20], p0 = [3, 40, 20], bounds = (0, np.inf))
params2, cov = curve_fit(lambda lam, B: B * light_inter(lam), x104[(y104 - ydark) >= 0][::50], (y104 - ydark)[(y104 - ydark) >= 0][::50], p0 = [40], bounds = (0, np.inf))
#params3, cov = curve_fit(sum_light_mie, x14[(y14 - ydark) >= 0][::50], (y14 - ydark)[(y14 - ydark) >= 0][::50], p0 = [3, 25, 20], bounds = (0, np.inf))


fig = plt.figure(figsize = (5.906, 4))
ax = fig.add_subplot(111)

ax.plot(x104, (y104 - ydark) , '.', label = '$\phi = 104$', alpha = 0.8, markersize = 0.8)
ax.plot(lam_plot, params1[-1] * light_inter(lam_plot), label = 'Lichtquelle')
ax.plot(lam_plot, params2[0] * light_inter(lam_plot), label = 'Nur Lichtquelle')
ax.plot(lam_plot, Q_sc(lam_plot, *params1[:-1]), 'k-', label = '$Q_{sc}$')
ax.plot(lam_plot, sum_light_mie(lam_plot, *params1), 'r-', label = 'Fit')
ax.set_xlim(lam_plot[0], lam_plot[-1])
ax.set_xlabel('Wellenlänge / nm')
ax.set_ylabel('Streu-Intensität / a.u.')
ax.legend()
fig.savefig(abs_path('results/fit_104.pdf'), bbox_inches = 'tight', pad_inches = 0)
