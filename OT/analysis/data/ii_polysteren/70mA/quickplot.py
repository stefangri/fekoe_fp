import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import os
name = os.getcwd().split('/')[-1]

x, y, sum = np.genfromtxt(f'{name}_ohne.dat', unpack = True)
xcal, ycal = np.genfromtxt(f'../../i_quarz/{name}/zcalibration.txt', unpack = True)
x /= xcal
y /= ycal
fs = 10e3

f, Pxx_den = signal.periodogram(x - np.mean(x), fs)
dt = 1 / fs
t = np.arange(0, len(x)*dt, step = dt)

fig = plt.figure(figsize = (5.7, 7.8))
ax1 = fig.add_subplot(211)
ax1.plot(f, Pxx_den)
ax1.set_xscale('log')
ax1.set_yscale('log')


f, Pxx_den = signal.periodogram(y - np.mean(y), fs)



ax2 = fig.add_subplot(212)
ax2.plot(f, Pxx_den)
ax2.set_xscale('log')
ax2.set_yscale('log')
plt.show()
fig.savefig('test.png')
