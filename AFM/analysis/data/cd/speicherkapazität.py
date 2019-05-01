#n Δx [10<sup>-6</sup>] Δy [10<sup>-6</sup>] φ [deg] R [10<sup>-6</sup>] Δz [10<sup>-9</sup>]
#1 1,533 -0,030 1,11 1,533 -27,68 #spurabstand
#2 0,020 -1,622 89,30 1,622 32,16 #max pitlänge
#3 0,010 -0,485 88,83 0,485 -5,76 #min pitlänge
#4 0,227 0,000 0,00 0,227 17,12 #spur breite

import math as m

spurabstand = 1.533
minpitlänge  = 0.485
print(minpitlänge / 2)
spurbreite  = 0.227

fläche =  	100 * 1e8

spurlänge = fläche / spurabstand
print(spurlänge * 1e-9)
num_bits =  spurlänge / minpitlänge
print(num_bits)
print('speicherkapzität: ', num_bits / (14 + 3) * 1e-9)
print(m.pi * (6**2 - 2.2**2))
