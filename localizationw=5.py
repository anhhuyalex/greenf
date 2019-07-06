import numpy as np
import time
import matplotlib.pyplot as plt
from collections import OrderedDict
import pickle
import datetime


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

def Hamiltonian(N, W = 1.0, disorder = True, e_n = 1.0, t_n = 1.0):
    """
        inputs
            N = side length of cubic lattice
            W = degree of disorder
        outputs
            H = Hamiltonian
            inv_mapping = coordinate of each state 
    """

    def gen_triples(N):
        i, j, k = (0, 0, 0)
        coords = []
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    coords.append((i,j,k))
        return coords

    coords = gen_triples(N)
    mapping = {i: coords[i] for i in range(len(coords))}
    inv_mapping = {coords[i]: i for i in range(len(coords))}

    sites = N**3
    H = np.zeros((sites, sites), dtype = complex)
    for i in range(sites):
        x, y, z = mapping[i]
        if disorder == True:
            H[i, i] = np.random.uniform(-W, W)
        else:
            H[i, i] = e_n
        if (z + 1 < N):
            H[i, i+1] = t_n
        if y + 1 < N:
            H[i, i+N] = t_n
        if x + 1 < N:
            H[i, i+N**2] = t_n
        if z > 0:
            H[i, i-1] = t_n
        if y > 0:
            H[i, i-N] = t_n
        if x > 0:
            H[i, i-N**2] = t_n
    return H, inv_mapping


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

    return new_map, new_inv_map

def _calculate(N, z_n, x_, y_, z_, Hamiltonian, site_coordinate, mapping, inv_mapping,
                t_n = 1.0, e_n = 1.0, k = 1, c = complex):

    
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

            # For Hamiltonians with order, substitute the correct amplitude
            site = site_coordinate[(x,y,z)]
            e_n = Hamiltonian[site, site]
            amplitude = t_n / (z_n - e_n)

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
    site = site_coordinate[(x_,y_,z_)]
    e_n = Hamiltonian[site, site]
    C[inv_mapping[sum_][(x_, y_, z_)]] = 1. / (z_n - e_n)
    
        



    # Calculate Multiplicative Factors 
    A = {}
    G = {}
    A[1] = b[1]
    tic = time.time()
    for k in range(2, sum_):
        # print("k", k)
        # print("a[k]", a[k].shape)
        # print("b[k]", b[k].shape)
        # print("A[k-1]", A[k - 1].shape)
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


def sign(control, target):
    if np.sign(control.real) == np.sign(target.real) and np.sign(control.imag) == np.sign(target.imag):
        return 1
    else:
        return -1

energy = 0 # np.random.uniform(-1, 1) 
imag = 0.1j
N = 20
W = 5
x_, y_, z_ = (19,19,18)
threshold = 0.01

print("Saving to file", 'greenfs-{}-{}-{}.pickle'.format(x_, y_, z_))



# true_value = true_value[sum([x_, y_, z_])][inv_mapping[sum([x_, y_, z_])][(0,0,0)]]


def drop_off(N, mapping, inv_mapping, x, y, z, x_, y_, z_, true_value, H, site_coordinate,
            threshold = threshold, energy=energy, imag=imag):
    xyzx_y_z_ = true_value[sum([x, y, z])][inv_mapping[sum([x, y, z])][(x, y, z)]]
    log_xyzx_y_z_ = np.log(np.absolute(xyzx_y_z_))
    omit = np.full(N**3, 1)
    ind = np.arange(N**3)
    batches = OrderedDict()
    batches[tuple(ind)] = 0 # keeps track of whether entire minibatch of states have been dropped

    continue_expand = True
    while continue_expand == True:
        new_batches = OrderedDict()
        print ([(len(batch), status) for batch, status in batches.items()])
        for batch, status in batches.items():
            if (status == 0) and len(batch) > 0: #if batch needs to be expanded
                
                minibatches = np.array_split(batch, 2) # divide batch of states to investigate into minibatches

                for minibatch in minibatches[::-1]:
                    
                    print ("len(minibatch)", len(minibatch))
                    omit[(minibatch)] = 0
                    _mapping, _inv_mapping = omit_sites(inv_mapping, omit)

                    if ((x_, y_, z_) not in _inv_mapping[sum([x_, y_, z_])]) or ((x, y, z) not in _inv_mapping[sum([x, y, z])]):
                        omit[(minibatch)] = 1
                        new_batches[tuple(minibatch)] = 0
                        continue
                    if (x_ + y_ + z_ == 3*N - 3) or (x_ + y_ + z_ == 0):
                        continue

                    predicted_value = (_calculate(N, energy+imag, x_, y_, z_, H, site_coordinate, _mapping, _inv_mapping))
                    predicted_value = predicted_value[sum([x, y, z])][inv_mapping[sum([x, y, z])][(x, y, z)]]
                    log_predicted_value = np.log(np.absolute(predicted_value)) * sign(predicted_value, xyzx_y_z_)

                    rel_err = np.absolute(( log_predicted_value - log_xyzx_y_z_  ) / log_xyzx_y_z_)
                    if rel_err > threshold:
                        omit[(minibatch)] = 1
                        new_batches[tuple(minibatch)] = 0
                    else:
                        new_batches[tuple(minibatch)] = 1
                    print ("true_value", xyzx_y_z_)
                    print ("predicted_value", predicted_value)

            else: # if batch doesn't affect accuracy, no need to expand
                new_batches[tuple(batch)] = 1
                omit[list(batch)] = 0

        batches = new_batches

        # Check if done with drop-off algorithm
        continue_expand = False
        for batch, status in batches.items():
            if (len(batch) > 1) and (status == 0):
                continue_expand = True
                break
    return predicted_value

for _ in range(15):
    mapping, inv_mapping = index_triple_maps(N)
    H, site_coordinate = Hamiltonian(N, W) #sample a Hamiltonian
    true_value = (_calculate(N, energy+imag, x_, y_, z_, H, site_coordinate, mapping, inv_mapping))
    greenf = OrderedDict()
    greenf["true_value"] = true_value

    R = [(mapping[i][0]) for i in true_value] # a certain element of distance R away
    # print(*R[-10:][::-1], sep='\n')
    predicted_value = OrderedDict()
    for x, y, z in R[-10:][::-1]:
        pred = drop_off(N, mapping, inv_mapping, x, y, z, x_, y_, z_, true_value, H, site_coordinate)
        predicted_value[sum([x, y, z])] = pred

    greenf["predicted_value"] = predicted_value
    now = datetime.datetime.now()
    root_dir = ""
    filename = "{}greenfs-w=5-{}.pickle".format(root_dir, now.isoformat())
    with open(filename, 'wb') as handle:
        pickle.dump(greenf, handle, protocol=pickle.HIGHEST_PROTOCOL)

# x, y, z = (0, 0, 1)

# true_value = true_value[sum([x, y, z])][inv_mapping[sum([x, y, z])][(x, y, z)]]
# 
# true_value = np.log(np.absolute(true_value)) 
# print("True value", true_value)

# print("batches", batches)
        
# plt.plot(range(8000), omit)
# plt.show()