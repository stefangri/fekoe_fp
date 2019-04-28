import numpy as np
from pint import UnitRegistry
u = UnitRegistry()
Q_ = u.Quantity

E = Q_(160, 'gigapascal')
w = Q_(50, 'micrometer')
l = Q_(450, 'micrometer')
t = Q_(2, 'micrometer')

k = E * w * t**3 / 4 / l**3
print(k.to('newton / meter'))
