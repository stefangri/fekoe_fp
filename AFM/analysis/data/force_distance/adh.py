import numpy as np
for name in ['edelstahl', 'teflon', 'dlc']:
    print(name)
    features = np.genfromtxt('features_force_distance_' + name + '.txt')
    pulloff = features[1]
    snapin = features[0]
    pulloff[0] *= (20 / 75)
    snapin[0] *= (20 / 75)

    k = 0.2 #N / m
    z = snapin[0] - pulloff[0] # mu

    F = k * z * 1e-6
    print(F)
