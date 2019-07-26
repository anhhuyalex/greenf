import numpy as np
import time
import matplotlib.pyplot as plt
from collections import OrderedDict
import pickle
import datetime
import localization

energy = 0 # np.random.uniform(-1, 1)
imag = 0.1j
N_small_lattice = 20
N_large_lattice = 30
W = 1.5
x_, y_, z_ = (1,2,2)
threshold = 0.01




print("Saving to file", 'greenfs-{}-{}-{}.pickle'.format(x_, y_, z_))

def extend_localization(energy = energy, imag = imag, N = N_small_lattice, W = W, x_ = x_, y_ = y_, z_ = z_,
                threshold = threshold):
    # for _ in range(15):
    print(N)
    mapping, inv_mapping = localization.index_triple_maps(N)
    H, site_coordinate = localization.Hamiltonian(N, W=0, disorder=False) #sample a Hamiltonian
    a, b = localization._build_ab_matrices(N, energy+imag, H, site_coordinate, mapping, inv_mapping)
    # true_value = (localization._calculate(N, energy+imag, x_, y_, z_, H, site_coordinate, mapping, inv_mapping))
    # for key in inv_mapping:
    #     print("Sum", key)
    #     for el in inv_mapping[key].keys():
    #         print("Elements: <", el, "|(1, 2, 2)>")
    # print(true_value)
    # greenf = OrderedDict()
    # greenf["true_value"] = true_value
    #
    # R = [(mapping[i][0]) for i in true_value] # a certain element of distance R away
    # # print(*R[-10:][::-1], sep='\n')
    # predicted_value = OrderedDict()
    # for x, y, z in R[-10:][::-1]:
    #     pred = drop_off(N, mapping, inv_mapping, x, y, z, x_, y_, z_, true_value, H, site_coordinate)
    #     predicted_value[sum([x, y, z])] = pred

    # greenf["predicted_value"] = predicted_value
    # now = datetime.datetime.now()
    # root_dir = ""
    # filename = "{}greenfs-w={}-{}.pickle".format(root_dir, str(W), now.isoformat())
    # with open(filename, 'wb') as handle:
    #     pickle.dump(greenf, handle, protocol=pickle.HIGHEST_PROTOCOL)

np.set_printoptions(precision=1)
extend_localization(N = 2)
extend_localization(N = 3)
