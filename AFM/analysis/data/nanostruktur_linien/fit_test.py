import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def lin(x, a, b):
    return a * x + b


x, y = np.genfromtxt('linien_profil.txt', unpack = True)

params, cov = curve_fit(lin, x, y)

plt.plot(x, y)
#plt.plot(x, lin(x, *params))
plt.show()
