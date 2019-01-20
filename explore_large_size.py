import numpy as np
import time
import matplotlib.pyplot as plt
import time
from collections import OrderedDict

energy = 0.8 # np.random.uniform(-1, 1) 
N=30
x_, y_, z_ = (5, 5, 5)
print("Saving to file", 'greenfs-{}-{}-{}.pickle'.format(x_, y_, z_))

mapping, inv_mapping = [], []

def generate_triples(N, total):
    # List containing generated triples
    gen = []

    
    if 0 <= total <= N - 1:
        # Generate first element
        i, j, k = (0, 0, total)
        while i <= total:
            # print ("now i is", i)
            while j <= total - i:
                gen.append ((i, j, k))
                j += 1
                k -= 1
            i += 1
            j = 0
            k = total - i
            
                # pass

    elif N <= total <= 2*N - 3:
        # Generate first element
        i, j, k = (0, total-N+1, N-1)
        while i <= N-1:
            while j <= min(N-1, total - i):
                gen.append ((i, j, k))
                j += 1
                k -= 1
            i += 1
            k = min(N-1, total - i)
            j = total - i - k


            
    elif 2*N - 2 <= total <= 3*N - 3:
        # Generate first element
        i, j, k = (total - 2*N + 2, N-1, N-1)
        while i <= N-1:
            while j <= min(N-1, total - i):
                gen.append ((i, j, k))
                j += 1
                k -= 1
            i += 1
            k = min(N-1, total - i)
            j = total - i - k
        
    else:
        raise ValueError("Not applicable k")

    return gen

def index_triple_maps(N):
    mapping = {}
    inv_mapping = {}
    for k in range(3*N - 2):
        triples = generate_triples(N, k)
        coords = OrderedDict()
        inv_coords = OrderedDict()
        for i in range(len(triples)):
            coords[i] = triples[i]
            inv_coords[triples[i]] = i
        mapping[k] = coords
        inv_mapping[k] = inv_coords
    return mapping, inv_mapping

def omit_sites(inv_mapping, omit):    
    site = 0
    new_inv_map = {}
    new_map = {}
    
    for sum_ in inv_mapping.keys():
        inv_sum_k_sites = OrderedDict()
        sum_k_sites = OrderedDict()
        site_index = 0
        for triple in inv_mapping[sum_].keys():
            if omit[site] == 1:
                inv_sum_k_sites[triple] = site_index
                sum_k_sites[site_index] = triple
                site_index += 1
            else:
                x, y, z = triple
                print ("I'm omitting", triple)
            site += 1
        new_inv_map[sum_] = inv_sum_k_sites
        new_map[sum_] = sum_k_sites

    return new_map, new_inv_map, x, y, z

def _calculate(z_n, x_, y_, z_, x_o, y_o, z_o, t_n = 1.0, e_n = 1.0, N = 20, k = 1, c = complex,
                           mapping = mapping, inv_mapping = inv_mapping):

    

    amplitude = t_n / (z_n - e_n)
    # amplitude = 1
    a = {}
    b = {}


    for k in range(3*N-2):
        alpha_num_entries = len(inv_mapping[k-1]) if k > 0 else 0
        my_entries = len(inv_mapping[k])
        beta_num_entries = len(inv_mapping[k+1]) if k < 3*N - 3 else 0
        a_k = np.zeros((my_entries, alpha_num_entries), dtype=c)
        b_k = np.zeros((my_entries, beta_num_entries), dtype=c)
        for index, triple in mapping[k].items():
            x, y, z = triple

            # Calculate b_k
            if k < 3*N - 3:
                if x < N - 1:
                    try:
                        b_k[index, inv_mapping[k + 1][x + 1, y, z]] = amplitude
                    except:
                        pass
                    
                if y < N - 1:
                    try:
                        b_k[index, inv_mapping[k + 1][x, y + 1, z]] = amplitude
                    except:
                        pass

                if z < N - 1:
                    try:
                        b_k[index, inv_mapping[k + 1][x, y, z + 1]] = amplitude
                    except:
                        pass
#             # Calculate a_k
            if 0 < k:
                if x > 0:
                    try:
                        a_k[index, inv_mapping[k - 1][x - 1, y, z]] = amplitude
                    except:
                        pass

                if y > 0:
                    try:
                        a_k[index, inv_mapping[k - 1][x, y - 1, z]] = amplitude
                    except:
                        pass

                if z > 0:
                    try:
                        a_k[index, inv_mapping[k - 1][x, y, z - 1]] = amplitude
                    except:
                        pass

        a[k] = a_k
        b[k] = b_k

    # Calculate C
    sum_ = sum([x_, y_, z_])
    
    C = np.zeros(shape = (len(inv_mapping[sum_]), 1), dtype = complex)
    C[inv_mapping[sum_][(x_, y_, z_)]] = 1. / (z_n - e_n)
    
        



    # Calculate Multiplicative Factors 
    A = {}
    G = {}
    A[1] = b[1]
    
    for k in range(2, sum_):
        A[k] = np.linalg.solve(np.eye(len(inv_mapping[k])) - a[k].dot(A[k - 1]), b[k])


    A[3*N - 4] = a[3*N - 4]
    for k in range(3*N - 5, sum_, -1):
        A[k] = np.linalg.solve(np.eye(len(inv_mapping[k])) - b[k].dot(A[k + 1]), a[k])

    
    # Calculate Green's Functions
    V = OrderedDict()
    k = sum_
    if 1 < k < 3*N-4:
        V[k] = np.linalg.solve(np.eye(len(inv_mapping[k])) - a[k].dot(A[k - 1]) - b[k].dot(A[k + 1]), 
                            C)
    elif k == 1:
        V[k] = np.linalg.solve(np.eye(len(inv_mapping[k])) - b[k].dot(A[k + 1]), 
                            C)
    elif k == 3*N - 4:
        V[k] = np.linalg.solve(np.eye(len(inv_mapping[k])) - a[k].dot(A[k - 1]), 
                            C)
            
    
    # if x_t + y_t + z_t > sum_:
    for k in range(sum_ + 1, 3*N - 3):
        V[k] = A[k].dot(V[k - 1])
    # elif x_t + y_t + z_t < sum_:
    for k in range(sum_ - 1, 0, -1):
        V[k] = A[k].dot(V[k + 1])

    # return V[x_t + y_t + z_t]
    return V


# print(inv_mapping[(1,4,2)])
# print (coords)


# Create training data

mapping, inv_mapping = index_triple_maps(N)
true_value = (_calculate(energy+0.1j, x_, y_, z_, 1, 1, 1, mapping = mapping, inv_mapping = inv_mapping))
pos = inv_mapping[sum([x_, y_, z_])][(x_, y_, z_)]
true_value = true_value[sum([x_, y_, z_])][pos]
greenfs = OrderedDict()  
greenfs["true_value"] = true_value
# X = []
# true_Y = OrderedDict()
pred_Y = OrderedDict()
# energy = np.random.uniform(-1, 1) 

print ("Position of element", pos)
for x_ in range(N):
    
    # x_, y_, z_ = np.random.randint(20, size=3)
    
    omit = np.full(N**3, 1)
    # for i in range(10):
        # omit[i] = 0
    omit[_] = 0
    _mapping, _inv_mapping, x_o, y_o, z_o = omit_sites(inv_mapping, omit)
    # print (_inv_mapping)
    # raise ValueError
    
    # Do not omit target and destination states
    if (x_ == x_o) & (y_ == y_o) & (z_ == z_o):
        continue
    if (x_ + y_ + z_ == 3*N - 3) or (x_ + y_ + z_ == 0):
        continue
    

    
    print (x_, y_, z_, x_o, y_o, z_o, energy)
    
    # true_value = (_calculate(energy+0.1j, x_, y_, z_, 1, 1, 1, mapping = mapping, inv_mapping = inv_mapping))

    # Remove omitted state from true value for error calculation
    # if x_o + y_o + z_o in true_value.keys():
    #     true_value[x_o + y_o + z_o] = np.delete(true_value[x_o + y_o + z_o], inv_mapping[x_o + y_o + z_o][(x_o, y_o, z_o)])
    #     true_value[x_o + y_o + z_o] = true_value[x_o + y_o + z_o][:,np.newaxis]

    
    predicted_value = (_calculate(energy+0.1j, x_, y_, z_, x_o, y_o, z_o, mapping = _mapping, inv_mapping = _inv_mapping))
    predicted_value = predicted_value[sum([x_, y_, z_])][_inv_mapping[sum([x_, y_, z_])][(x_, y_, z_)]]
    
    # Calculate error 
    # for k in true_value.keys():
        # relative error = modulus of error / modulus of true value
        # error = np.absolute((predicted_value[k] - true_value[k]) / true_value[k] * 100)
        # print ("I'm printing key ", k)
#         print ("I'm printing true_value")
#         print (true_value[k])
        # print ("I'm printing error")
#         print (predicted_value[k])
        # print (np.absolute(true_value[k] - predicted_value[k]) / np.absolute(true_value[k]))
        # for index, triple in enumerate(inv_mapping[k].keys()):
        #     (x, y, z) = triple
        #     x, y, z, x__, y__, z__, x_o_, y_o_, z_o_ = (np.array([x, y, z, x_, y_, z_, x_o, y_o, z_o], dtype=float) / 10.0) - 1
            # x, y, z, x__, y__, z__ = (np.array([x, y, z, x_, y_, z_], dtype=int)) 
        #     X.append([x, y, z, x__, y__, z__, energy])
        # print(true_value[k])
        # print(predicted_value[k])
        # true_Y.extend(true_value[k])
        # pred_Y.extend(predicted_value[k])
    # true_Y[_] = true_value
    pred_Y[_] = predicted_value
    print ("true_value", true_value)
    print ("predicted_value", predicted_value)
    # print (len(X))

greenfs["predicted_value"] = pred_Y
# print (true_Y)
# print (true_Y.keys())
# print (true_Y[1000].keys())
# print (pred_Y)
import pickle
with open('greenfs-{}-{}-{}.pickle'.format(x_, y_, z_), 'wb') as handle:
    pickle.dump(greenfs, handle, protocol=pickle.HIGHEST_PROTOCOL)
# with open('pred_second.pickle', 'wb') as handle:
#     pickle.dump(pred_Y, handle, protocol=pickle.HIGHEST_PROTOCOL)

# true_Y = np.array(true_Y)
# pred_Y = np.array(pred_Y)
# X = np.array(X)
# np.save("predictors2" + str(energy), X)
# np.save("true", true_Y)
# np.save("pred", pred_Y)