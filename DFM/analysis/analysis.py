import numpy as np
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as stds
from uncertainties import correlated_values
import math
from scipy.optimize import curve_fit
from pint import UnitRegistry
import matplotlib.gridspec as gridspec
import latex as l
u = UnitRegistry()
Q_ = u.Quantity
import os

def abs_path(filename):
    return os.path.join(os.path.dirname(__file__), filename)

def fwhm(lam, Q_sc):
    max_wl = lam[np.argmax(Q_sc)]
    half_above = lam[np.argmax(Q_sc):][np.argmin(abs(Q_sc[np.argmax(Q_sc):] - max(Q_sc)/2))]
    half_below = lam[:np.argmax(Q_sc)][np.argmin(abs(Q_sc[:np.argmax(Q_sc)] - max(Q_sc)/2))]
    return {'max_I': max(Q_sc), 'max_wl': max_wl, 'lam+': half_above, 'lam-': half_below, 'fwhm': half_above - half_below}

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




#l, ref_n, k = np.genfromtxt(abs_path('data/refrective_index.txt'), unpack = True) #data for refrective index
k, l, ref_n = np.genfromtxt(abs_path('data/refrective_index2.txt'), unpack = True) #data for refrective index
l *= 1e3 #convert lambda to nanometer


n_inter = interpolate.interp1d(l, ref_n, kind = 'quadratic') #interpolate real part of refrective index
k_inter = interpolate.interp1d(l, k, kind = 'quadratic') #interpolate imaginary part of refrective index

lam_plot = np.linspace(400, 800, 1000) #wavelength array for plots

fig = plt.figure(figsize = (5.906, 2.5))
ax1 = fig.add_subplot(121)
ax1.text(0.01 , 0.96, r'(a)', horizontalalignment='left', verticalalignment='top', transform = ax1.transAxes, color = 'black')
ax1.plot(l[(l >=400) & (l <= 800)], ref_n[(l >=400) & (l <= 800)], '.', label = 'Messwerte')
ax1.plot(lam_plot, n_inter(lam_plot), label = 'Interpolation')
ax1.set_xlabel('Wellenlänge')
ax1.set_ylabel('$n$')
ax1.set_xlim(lam_plot[0], lam_plot[-1])
ax1.legend()

ax2 = fig.add_subplot(122)
ax2.text(0.01 , 0.96, r'(b)', horizontalalignment='left', verticalalignment='top', transform = ax2.transAxes, color = 'black')
ax2.plot(l[(l >=400) & (l <= 800)], k[(l >=400) & (l <= 800)], '.')
ax2.plot(lam_plot, k_inter(lam_plot))
ax2.set_xlabel('Wellenlänge')
ax2.set_ylabel('$k$')
ax2.set_xlim(lam_plot[0], lam_plot[-1])

fig.tight_layout()
fig.savefig(abs_path('results/refrective_index.png'), bbox_inches = 'tight', pad_inches = 0)



x104, y104 = np.genfromtxt(abs_path('data/sag5b_0001.txt'), unpack = True) #load spectrum data
x14, y14 = np.genfromtxt(abs_path('data/sag6_0001.txt'), unpack = True)
xdark, ydark = np.genfromtxt(abs_path('data/sag7_0001.txt'), unpack = True)
pic104 = mpimg.imread(abs_path('data/pics/au_röhren_spektrum_02_sag5.png'))
pic14  = mpimg.imread(abs_path('data/pics/au_röhren_spektrum_02_sag6.png'))


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



params104, cov   = curve_fit(sum_light_mie, x104[::20], (y104 - ydark)[::20], p0 = [3, 40, 20], bounds = (0, np.inf))
paramslight, cov = curve_fit(lambda lam, B: B * light_inter(lam), x104[::20], (y104 - ydark)[::20], p0 = [40], bounds = (0, np.inf))
print(params104[1])
print(paramslight)

fig = plt.figure(figsize = (5.906, 4 * 1.5))
spec = gridspec.GridSpec(ncols=2, nrows=3, figure=fig)
axpic104 = fig.add_subplot(spec[0, 0])
axpic104.text(0.01 , 0.96, r'(a) $\phi = \SI{104}{\degree}$', horizontalalignment='left', verticalalignment='top', transform = axpic104.transAxes, color = 'white')
axpic104.imshow(pic104)
axpic104.set_yticks([])
axpic104.set_xticks([])

axpic14 = fig.add_subplot(spec[0, 1])
axpic14.text(0.01 , 0.96, r'(b) $\phi = \SI{14}{\degree}$', horizontalalignment='left', verticalalignment='top', transform = axpic14.transAxes, color = 'white')
axpic14.imshow(pic14)
axpic14.set_yticks([])
axpic14.set_xticks([])

ax = fig.add_subplot(spec[1:, 0])

ax.text(0.01 , 0.96, r'(c)', horizontalalignment='left', verticalalignment='top', transform = ax.transAxes, color = 'black')
ax.plot(x104, (y104 - ydark) , '.', label = r'$\phi = \SI{104}{\degree}$', alpha = 0.8, markersize = 0.8)
ax.plot(lam_plot, Q_sc(lam_plot, *params104[:-1]), 'k-', label = '$Q_{sc}$')
ax.plot(lam_plot, params104[-1] * light_inter(lam_plot), label = 'Lichtquelle')
ax.plot(lam_plot, sum_light_mie(lam_plot, *params104), 'r-', label = 'Summe')
ax.plot(lam_plot, paramslight[0] * light_inter(lam_plot), label = 'Nur Lichtquelle')
ax.set_xlim(lam_plot[0], lam_plot[-1])
ax.set_xlabel('Wellenlänge / nm')
ax.set_ylabel('Streu-Intensität / a.u.')
ax.legend()
#fig.savefig(abs_path('results/fit_104.png'), bbox_inches = 'tight', pad_inches = 0)



params14, cov    = curve_fit(sum_light_mie, x14[::20], (y14 - ydark)[::20], p0 = [3, 80, 30], bounds = (0, np.inf))
paramslight, cov = curve_fit(lambda lam, B: B * light_inter(lam), x14[::20], (y14 - ydark)[::20], p0 = [40], bounds = (0, np.inf))
print(params14[1])
print(paramslight)

#fig = plt.figure(figsize = (5.906, 4))
ax2 = fig.add_subplot(spec[1:, 1], sharey = ax)
ax2.text(0.01 , 0.96, r'(d)', horizontalalignment='left', verticalalignment='top', transform = ax2.transAxes, color = 'black')
ax2.plot(x14, (y14 - ydark) , '.', label = r'$\phi = \SI{14}{\degree}$', alpha = 0.8, markersize = 0.8)
ax2.plot(lam_plot, params14[-1] * light_inter(lam_plot))
ax2.plot(lam_plot, paramslight[0] * light_inter(lam_plot))
ax2.plot(lam_plot, Q_sc(lam_plot, *params14[:-1]), 'k-')
ax2.plot(lam_plot, sum_light_mie(lam_plot, *params14), 'r-')
ax2.set_xlim(lam_plot[0], lam_plot[-1])
ax2.set_xlabel('Wellenlänge / nm')
#ax2.set_ylabel('Streu-Intensität / a.u.')
ax2.legend(loc = 'upper right')
fig.savefig(abs_path('results/fit_14_104.png'), bbox_inches = 'tight', pad_inches = 0)




fig = plt.figure(figsize = (5.906, 2.5))
ax = fig.add_subplot(111)

peak_info_14 = fwhm(lam_plot, Q_sc(lam_plot, *params14[:-1]))
peak_info_104 = fwhm(lam_plot, Q_sc(lam_plot, *params104[:-1]))

print(peak_info_14)
print(peak_info_104)


ax.plot(lam_plot, Q_sc(lam_plot, *params14[:-1]), label = r'$Q_{sc}, \phi = \SI{14}{\degree}$', color = 'b')
ax.plot([peak_info_14['lam-'], peak_info_14['lam+']], [peak_info_14['max_I']/2, peak_info_14['max_I']/2], 'b--')
ax.plot(lam_plot, Q_sc(lam_plot, *params104[:-1]), label = r'$Q_{sc}, \phi = \SI{104}{\degree}$', color = 'r')
ax.plot([peak_info_104['lam-'], peak_info_104['lam+']], [peak_info_104['max_I']/2, peak_info_104['max_I']/2], 'r--')
ax.set_xlim(lam_plot[0], lam_plot[-1])
ax.set_xlabel('Wellenlänge / nm')
ax.set_ylabel('Streu-Intensität / a.u.')
ax.legend()
fig.savefig(abs_path('results/compare_14_104.png'), bbox_inches = 'tight', pad_inches = 0)




xspheres, yspheres = np.genfromtxt(abs_path('data/sag8_0001.txt'), unpack = True) #load spectrum data
xdark, ydark = np.genfromtxt(abs_path('data/sag9_0001.txt'), unpack = True)


#paramsspheres, cov   = curve_fit(sum_light_mie, xspheres[::20], (yspheres)[::20], p0 = [15, 100, 100], bounds = (0, np.inf))
#paramslight, cov = curve_fit(lambda lam, B: B * light_inter(lam), xspheres[::20], (yspheres)[::20], p0 = [40], bounds = (0, np.inf))


fig = plt.figure(figsize = (5.906, 2.5))

ax = fig.add_subplot(131)
ax.plot(xspheres, yspheres , '.', alpha = 0.8, markersize = 0.8)
ax.set_xlim(lam_plot[0], lam_plot[-1])
ax.text(0.01 , 0.96, r'(a) Nanosphären', horizontalalignment='left', verticalalignment='top', transform = ax.transAxes, color = 'black')
ax.set_ylabel('Streu-Intensität / a.u.')
#ax.legend()


ax2 = fig.add_subplot(132, sharey = ax)
ax2.plot(xdark, ydark , '.', alpha = 0.8, markersize = 0.8)
ax2.set_xlim(lam_plot[0], lam_plot[-1])
ax2.text(0.01 , 0.96, r'(b) Dunkel', horizontalalignment='left', verticalalignment='top', transform = ax2.transAxes, color = 'black')
ax2.set_xlabel('Wellenlänge / nm')
#ax2.legend()



ax3 = fig.add_subplot(133)
ax3.plot(xspheres, yspheres - ydark, '.', alpha = 0.8, markersize = 0.8)
ax3.set_xlim(lam_plot[0], lam_plot[-1])
ax3.text(0.01 , 0.96, r'(c) Differenz', horizontalalignment='left', verticalalignment='top', transform = ax3.transAxes, color = 'black')#, bbox={'facecolor':'white', 'alpha':0.5, 'pad':0, 'edgecolor': None}
fig.tight_layout()
fig.savefig(abs_path('results/fit_spheres.png'), bbox_inches = 'tight', pad_inches = 0)




#Plot with CMOS pictures
fig = plt.figure(figsize = (5.906, 4))
pics = ['au_röhren_04.png', 'au_röhren_06.png', 'au_röhren_07.png', 'smiley2.png']
for i in range(4):
    pic = mpimg.imread(abs_path('data/pics/' + pics[i]))
    ax = fig.add_subplot(221 + i)
    ax.text(0.01 , 0.96, f'(' + ['a', 'b', 'c', 'd'][i] + ')', horizontalalignment='left', verticalalignment='top', transform = ax.transAxes, color = 'white')
    ax.imshow(pic)
    ax.set_yticks([])
    ax.set_xticks([])
fig.savefig(abs_path('results/pics_röhren.png'), bbox_inches = 'tight', pad_inches = 0)


fig = plt.figure(figsize = (5.906, 4))
pics = ['au_sphären_01.png', 'au_sphären_02.png']
for i in range(len(pics)):
    pic = mpimg.imread(abs_path('data/pics/' + pics[i]))
    ax = fig.add_subplot(121 + i)
    ax.text(0.01 , 0.96, f'(' + ['a', 'b'][i] + ')', horizontalalignment='left', verticalalignment='top', transform = ax.transAxes, color = 'white')
    ax.imshow(pic)
    ax.set_yticks([])
    ax.set_xticks([])
fig.savefig(abs_path('results/pics_spheres.png'), bbox_inches = 'tight', pad_inches = 0)
