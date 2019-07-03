import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import glob
from scipy.optimize import curve_fit

def power_density(f, A, f0):
    return A / (f**2 + f0**2)


dat_files = np.sort(glob.glob("*/*x.dat"))

for file in dat_files:
    x, y, sum = np.genfromtxt(file, unpack = True)
    foldername = file.split('/')[0]
    xcal, ycal = np.genfromtxt(f'{foldername}/zcalibration.txt', unpack = True)
    x /= xcal
    y /= ycal
    fs = 10e3
    f, Pxx_den = signal.periodogram(x-np.mean(x), fs)
    dt = 1 / fs
    t = np.arange(0, len(x)*dt, step = dt)
    plt.clf()
    plt.plot(t, x)
    plt.xlabel('t')
    plt.ylabel('x / $\mu$m')
    plt.savefig(f'{foldername}/xplot_{foldername}.png')
    plt.clf()

    plt.plot(f, Pxx_den)
    Pxx_den[f < 10] = Pxx_den[f>10][0]
    x, y, z, g = np.genfromtxt(f'{foldername}/force_calibration_{foldername}.FDdat', unpack = True)

    plt.plot(x, z)#[(x > 10) & (x < 3000)]

    params, cov = curve_fit(power_density, f[(f > 10) & (f < 800)], Pxx_den[(f > 10) & (f < 800)], p0 = [1e-4, 2])#,
    plt.axvline(x = params[1])
    print(foldername, params)

    fplot = np.logspace(-1, 5)

    params2, cov2 = curve_fit(power_density, x[(x > 10) & (x < 800)], z[(x > 10) & (x < 800)], p0 = [1e-8, 10])
    print(params2)


    plt.plot(fplot, power_density(fplot, *params2), 'b-')
    plt.plot(fplot, power_density(fplot, *params), 'r-')
    plt.xlabel('Frequenz')
    plt.ylabel('Power density')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([1e-17, 1e-3])
    plt.savefig(f'{foldername}/xpowerdensity_{foldername}.png')


dat_files = glob.glob("*/*y.dat")

for file in dat_files:
    x, y, sum = np.genfromtxt(file, unpack = True)
    dt = 1 / 10000
    foldername = file.split('/')[0]
    xcal, ycal = np.genfromtxt(f'{foldername}/zcalibration.txt', unpack = True)
    x /= xcal
    y /= ycal
    t = np.arange(0, len(y)*dt, step = dt)
    plt.clf()
    plt.plot(t, y)
    plt.xlabel('t')
    plt.ylabel('y /$\mu$m')
    plt.savefig(f'{foldername}/yplot_{foldername}.png')


dat_files = glob.glob("*/*ohne.dat")

for file in dat_files:
    x, y, sum = np.genfromtxt(file, unpack = True)
    foldername = file.split('/')[0]
    xcal, ycal = np.genfromtxt(f'{foldername}/zcalibration.txt', unpack = True)
    x /= xcal
    y /= ycal
    dt = 1 / 10000
    t = np.arange(0, len(x)*dt, step = dt)
    plt.clf()
    plt.plot(t, x - np.mean(x), label = 'x')
    plt.plot(t, y - np.mean(y), label = 'y')
    plt.xlabel('t')
    plt.ylabel('$\mu$m')
    plt.savefig(f'{foldername}/ohneplot_{foldername}.png')
