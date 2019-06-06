import numpy as np
import scipy.constants as const
# import pint
# ureg = pint.UnitRegistry()


me = 0.13*const.m_e
Eg = 1.84*const.e
print(Eg)
# mh = -1*const.m_e  # parallel
mh = -0.45*const.m_e  # senkrecht

# e_r = 9.29  # parallel
e_r = 9.15  # senkrecht

Epl1 = const.h*(1/(520*10**(-9)))*const.c
Epl2 = const.h*(1/(580*10**(-9)))*const.c
Epl3 = const.h*(1/(645*10**(-9)))*const.c
# print(Eg-Epl1)
# print(Epl1)
# print(Epl3)
mu = me*mh/(mh+me)

print(const.m_e)
c1 = ((const.hbar * const.pi)**2) * (1/mh + 1/me)
b1 = Eg - Epl1 - ((mu * const.e**4))/ 2 * (4 * const.pi * const.epsilon_0 * e_r * const.hbar)**2


q1 = c1 / b1
p1 = (1.786*const.e**2) / (8 * const.pi * e_r * const.epsilon_0 * (Eg - Epl1 - (mu * const.e**4) / (2*(4*const.pi*const.epsilon_0 * e_r * const.hbar)**2)))

a1 = p1/2 + np.sqrt(p1**2 * 1/4 - q1)
print(a1)

c2 = ((const.hbar * const.pi)**2) * (1/mh + 1/me)
b2 = Eg - Epl2 - ((mu * const.e**4)) / 2 * (4 * const.pi * const.epsilon_0 * e_r * const.hbar)**2


q2 = c2 / b2
p2 = (1.786*const.e**2) / (8 * const.pi * e_r * const.epsilon_0 * (Eg - Epl2 - (mu * const.e**4) / (2*(4*const.pi*const.epsilon_0 * e_r * const.hbar)**2)))

a2 = p2/2 + np.sqrt(p2**2 * 1/4 - q2)
print(a2)

c3 = ((const.hbar * const.pi)**2) * (1/mh + 1/me)
b3 = Eg - Epl3 - ((mu * const.e**4)) / 2 * (4 * const.pi * const.epsilon_0 * e_r * const.hbar)**2


q3 = c3 / b3
p3 = (1.786*const.e**2) / (8 * const.pi * e_r * const.epsilon_0 * (Eg - Epl3 - (mu * const.e**4) / (2*(4*const.pi*const.epsilon_0 * e_r * const.hbar)**2)))

a3 = p3/2 + np.sqrt(p3**2 * 1/4 - q3)
print(a3)


# B1 = ((const.hbar * const.pi)**2)/(2*(Eg-Epl1)) * (1/mh + 1/me) + (mu * const.e**4) / (2*(4*const.pi*const.epsilon_0 * e_r * const.hbar)**2*(Eg-Epl1))
# B2 = ((const.hbar * const.pi)**2)/(2*(Eg-Epl2)) * (1/mh + 1/me) + (mu * const.e**4) / (2*(4*const.pi*const.epsilon_0 * e_r * const.hbar)**2*(Eg-Epl2))
# B3 = ((const.hbar * const.pi)**2)/(2*(Eg-Epl3)) * (1/mh + 1/me) + (mu * const.e**4) / (2*(4*const.pi*const.epsilon_0 * e_r * const.hbar)**2*(Eg-Epl3))
# a1 = (1.786*const.e**2) / (8 * const.pi * e_r * const.epsilon_0 * (Eg - Epl1)) + np.sqrt(((1.786 * const.e**2) / (8 * const.pi * const.epsilon_0*e_r*(Eg-Epl1)))**2 * (1/4)-B1)
# a2 = (1.786*const.e**2) / (8 * const.pi * e_r * const.epsilon_0 * (Eg - Epl2)) + np.sqrt(((1.786 * const.e**2) / (8 * const.pi * const.epsilon_0*e_r*(Eg-Epl2)))**2 * (1/4)-B2)
# a3 = (1.786*const.e**2) / (8 * const.pi * e_r * const.epsilon_0 * (Eg - Epl3)) + np.sqrt(((1.786 * const.e**2) / (8 * const.pi * const.epsilon_0*e_r*(Eg-Epl3)))**2 * (1/4)-B3)

# print()
# print(const.epsilon_0)
# print()
# # print('effektive massen und epsilon_r parallel:')
# print('effektive massen und epsilon_r senkrecht:')
# print()
# print('lambda = 520nm : r = ', a1)
# print('lambda = 580nm : r = ', a2)
# print('lambda = 645nm : r = ', a3)
