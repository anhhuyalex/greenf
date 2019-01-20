import numpy as np
import time
import matplotlib.pyplot as plt
from collections import OrderedDict

energy = 0.8 # np.random.uniform(-1, 1) 
imag = 0.1j
N=30
x_, y_, z_ = (11,16,13)
print("Saving to file", 'greenfs-{}-{}-{}.pickle'.format(x_, y_, z_))

mapping, inv_mapping = [], []

import sys

# class Logger(object):
#     def __init__(self):
#         self.terminal = sys.stdout
#         self.log = open("log.dat", "a")

#     def write(self, message):
#         self.terminal.write(message)
#         self.log.write(message)  

# sys.stdout = Logger()

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
    print ("Omitting", (omit).shape[0] - np.sum(omit), "sites")    
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
                # print ("I'm omitting", triple)
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

    # tic = time.time()
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
    # toc = time.time()
    # print ("Time building a, b matrices", toc - tic)

    # Calculate C
    sum_ = sum([x_, y_, z_])
    
    C = np.zeros(shape = (len(inv_mapping[sum_]), 1), dtype = complex)
    C[inv_mapping[sum_][(x_, y_, z_)]] = 1. / (z_n - e_n)
    
        



    # Calculate Multiplicative Factors 
    A = {}
    G = {}
    A[1] = b[1]
    tic = time.time()
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
    toc = time.time()
    print ("Time doing lin alg", toc - tic)
    return V


mapping, inv_mapping = index_triple_maps(N)
true_value = (_calculate(energy+imag, x_, y_, z_, 1, 1, 1, mapping = mapping, inv_mapping = inv_mapping))
true_value = true_value[sum([x_, y_, z_])][inv_mapping[sum([x_, y_, z_])][(x_, y_, z_)]]
greenfs = OrderedDict()  
greenfs["true_value"] = true_value
threshold = 0.01

from collections import OrderedDict

omit = np.full(N**3, 1)

ind = np.arange(N**3)
batches = OrderedDict()
batches[tuple(ind)] = 0

minibatch_size = 600
tic = time.time()
while minibatch_size >= 1:
    # if minibatch_size < 10:
    #     threshold = 0.02
    # elif minibatch_size < 3:
    #     threshold = 0.01
    new_batches = OrderedDict()
    print ([(len(batch), status) for batch, status in batches.items()])
    print ("minibatch_size", minibatch_size)
    for batch, status in batches.items():
        if status == 0:
            dropping = len(batch) // minibatch_size
            
            minibatches = np.array_split(batch, dropping)
            print ("dropping", dropping)
    #         print(np.array_split(batch, dropping))

            for minibatch in minibatches[::-1]:
                
                print ("len(minibatch)", len(minibatch))
                omit[(minibatch)] = 0
                _mapping, _inv_mapping, x_o, y_o, z_o = omit_sites(inv_mapping, omit)

                if ((x_, y_, z_) not in _inv_mapping[sum([x_ + y_ + z_])]):
                    omit[(minibatch)] = 1
                    new_batches[tuple(minibatch)] = 0
                    continue
                if (x_ + y_ + z_ == 3*N - 3) or (x_ + y_ + z_ == 0):
                    continue

                predicted_value = (_calculate(energy+imag, x_, y_, z_, x_o, y_o, z_o, mapping = _mapping, inv_mapping = _inv_mapping))
                predicted_value = predicted_value[sum([x_, y_, z_])][_inv_mapping[sum([x_, y_, z_])][(x_, y_, z_)]]

                rel_err = np.absolute(predicted_value - true_value) / np.absolute(true_value)
                if rel_err > threshold:
                    omit[(minibatch)] = 1
                    print ("Error", rel_err, "over threshold!", threshold)
                    new_batches[tuple(minibatch)] = 0
                else:
                    new_batches[tuple(minibatch)] = 1
                print ("true_value", true_value)
                print ("predicted_value", predicted_value)

        elif status == 1:
            new_batches[tuple(batch)] = 1
            omit[list(batch)] = 0
    
    minibatch_size = minibatch_size / 4
    if 1 < minibatch_size < 2:
        minibatch_size = 1
    batches = new_batches
toc = time.time()
print("Calculation took", toc-tic)
        
        