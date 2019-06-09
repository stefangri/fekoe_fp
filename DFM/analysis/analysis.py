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
r = l.Latexdocument('results.tex')
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
from scipy.special import riccati_yn
from scipy.special import riccati_jn
from scipy.optimize import curve_fit
from scipy import interpolate

xlight, ylight = np.genfromtxt('light2.txt', unpack = True)
xlight = xlight[1:]
ylight = ylight[1:] / max(ylight[1:])
light_inter = interpolate.interp1d(xlight, ylight, kind = 'linear')




l, ref_n, k = np.genfromtxt('refrective_index.txt', unpack = True)
l *= 1e3

lam_plot = np.linspace(400, 800, 1000)
#n_poly = np.poly1d(np.polyfit(l, ref_n, 20))

n_inter = interpolate.interp1d(l, ref_n, kind = 'quadratic')
#plt.plot(lam_plot, n_inter(lam_plot))
#plt.plot(l, ref_n, 'ro')
#plt.show()

k_poly = np.poly1d(np.polyfit(l, k, 10)) #[(l>=400) & (l <= 800)]
k_inter = interpolate.interp1d(l, k, kind = 'quadratic')
#plt.plot(lam_plot, k_inter(lam_plot))
#plt.plot(l, k, 'ro')
#plt.show()

x104, y104 = np.genfromtxt('sag5b_0001.txt', unpack = True)
x14, y14 = np.genfromtxt('sag6_0001.txt', unpack = True)
xdark, ydark = np.genfromtxt('sag7_0001.txt', unpack = True)



N_medium = 1.5
#N_partikel = (0.2 + 3j)


def x(lam, size):
    return  np.pi / lam  * N_medium * size



def a_n(x, m, n):
    return np.array([(m * riccati_jn(n, m * x)[0][-1] * riccati_jn(n, x)[1][-1] - riccati_jn(n, x)[0][-1] * riccati_jn(n, m * x)[1][-1]) / (m * riccati_jn(n, m * x)[0][-1] * riccati_yn(n, x)[1][-1] - riccati_jn(n, m * x)[1][-1] * riccati_yn(n, x)[0][-1]) for x, m in zip(x, m)])


def b_n(x, m, n):
    return np.array([(riccati_jn(n, m * x)[0][-1] * riccati_jn(n, x)[1][-1] - m * riccati_jn(n, x)[0][-1] * riccati_jn(n, m * x)[1][-1]) / (riccati_jn(n, m * x)[0][-1] * riccati_yn(n, x)[1][-1] - m * riccati_jn(n, m * x)[1][-1] * riccati_yn(n, x)[0][-1]) for x, m in zip(x, m)])





#def m(lam)






def Q_sc(lam, A, size):
    #lam_mask = np.array([np.argmin(np.abs(l - lam)) for lam in lam])
    #N_partikel = np.array([ref_n[index]*(1 + 0j) + k[index]*(0 + 1j) for index in lam_mask])
    #N_partikel = n_poly(lam) * (1 + 0j) + k_poly(lam)*(0 + 1j)
    #N_partikel = np.interp(lam, l, ref_n) * (1 + 0j) + np.interp(lam, l, k)*(0 + 1j)
    n = 1
    N_partikel = n_inter(lam) * (1 + 0j) + k_inter(lam)*(0 + 1j)
    m = N_partikel / N_medium
    x_value = x(lam, size)
    a = np.array([a_n(x_value, m, n) for n in range(1, n + 1)])
    b = np.array([b_n(x_value, m, n) for n in range(1, n + 1)])
    return A * 2 / x_value**2 * np.sum( (2 * n + 1) * np.real(a * np.conj(a) + b * np.conj(b)), axis = 0)

#def Q_ext(lam, A, size, n = 1):
#    N_partikel = n_poly(lam)*(1 + 0j) + k_poly(lam)*(0 + 1j)
#    m = N_partikel / N_medium
#    x_value = x(lam, size)
#    a = np.array([a_n(x_value, m, n) for n in range(1, n + 1)])
#    b = np.array([b_n(x_value, m, n) for n in range(1, n + 1)])
#    return A * 2 / x_value**2  * np.sum( (2 * n + 1) * np.real(a + b), axis = 0)

def sum_light_mie(lam, A, size, B):
    return Q_sc(lam, A, size) + B * light_inter(lam)

#plt.plot(lam_plot, Q_sc(lam_plot, 1,60))
##plt.plot(lam_plot, -Q_ext(lam_plot, 1, 300))
#plt.show()

params1, cov = curve_fit(sum_light_mie, x104[(y104 - ydark) >= 0][::50], (y104 - ydark)[(y104 - ydark) >= 0][::50], p0 = [30, 300, 40], bounds = (0, np.inf))
params2, cov = curve_fit(lambda lam, B: B * light_inter(lam), x104[(y104 - ydark) >= 0][::50], (y104 - ydark)[(y104 - ydark) >= 0][::50], p0 = [40], bounds = (0, np.inf))
print(params1)




plt.plot(x104, (y104 - ydark) , ',', label = '$\phi = 104$', alpha = 1)
plt.plot(lam_plot, params1[-1] * light_inter(lam_plot), label = 'Lichtquelle')
plt.plot(lam_plot, params2[0] * light_inter(lam_plot), label = 'Nur Lichtquelle')
plt.plot(lam_plot, Q_sc(lam_plot, *params1[:-1]), 'k-', label = 'Nanoröhren')
plt.plot()
plt.plot(lam_plot, sum_light_mie(lam_plot, *params1), 'r-', label = 'Fit')
plt.xlim(lam_plot[0], lam_plot[-1])
plt.xlabel('Wellenlänge / nm')
plt.ylabel('Streu-Intensität / a.u.')
plt.legend()
#plt.show()
plt.savefig('test.png')
