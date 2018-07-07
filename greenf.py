import numpy as np
import time
import matplotlib.pyplot as plt
import time

tic = time.time()
def num_entries(k, N):
    if 0 <= k <= N - 1:
        return int(1./2 * (k+1) * (k+2))
    elif N <= k <= 2*N - 3:
        return int(1./2 * (-2*k**2 + 6*k*N - 6*k - 3*N**2 +9*N - 4))
    elif 2*N - 2 <= k <= 3*N - 3:
        return int(1./2 * (3*N - 2 - k) * (3*N - 1 - k))
    else:
        return 0


def generate_triples(N, total):
    # List containing generated triples
    gen = []

    # Generate first element
    if 0 <= total <= N - 1:
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
        mapping[k] = {i: triples[i] for i in range(len(triples))}
        inv_mapping[k] = {triples[i]: i for i in range(len(triples))}
    return mapping, inv_mapping



N = 1000

mapping, inv_mapping = index_triple_maps(N)
def calculate(z_n, t_n = 1.0, e_n = 1.0, N = 20, k = 1, c = complex,
                           mapping = mapping, inv_mapping = inv_mapping):

    

    amplitude = t_n / (z_n - e_n)
    # amplitude = 1
    a = {}
    b = {}


    tic = time.time()
    for k in range(3*N-2):
        alpha_num_entries = num_entries(k - 1, N)
        my_entries = num_entries(k, N)
        beta_num_entries = num_entries(k + 1, N)
        a_k = np.zeros((my_entries, alpha_num_entries), dtype=c)
        b_k = np.zeros((my_entries, beta_num_entries), dtype=c)
        for index, triple in mapping[k].items():
            x, y, z = triple

            # Calculate b_k
            if k < 3*N - 3:
                if x < N - 1:
                    b_k[index, inv_mapping[k + 1][x + 1, y, z]] = amplitude

                if y < N - 1:
                    b_k[index, inv_mapping[k + 1][x, y + 1, z]] = amplitude

                if z < N - 1:
                    b_k[index, inv_mapping[k + 1][x, y, z + 1]] = amplitude

            # Calculate a_k
            if 0 < k:
                if x > 0:
                    a_k[index, inv_mapping[k - 1][x - 1, y, z]] = amplitude

                if y > 0:
                    a_k[index, inv_mapping[k - 1][x, y - 1, z]] = amplitude

                if z > 0:
                    a_k[index, inv_mapping[k - 1][x, y, z - 1]] = amplitude

        a[k] = a_k
        b[k] = b_k

    # Calculate C
    x_, y_, z_ = (500,516,601)
    sum_ = sum([x_, y_, z_])
    C = np.zeros(shape = (num_entries(sum_, N), 1), dtype = complex)
    C[inv_mapping[sum_][(x_, y_, z_)]] = 1. / (z_n - e_n)



    # Calculate Multiplicative Factors 
    A = {}
    G = {}
    A[1] = b[1]
    # G[1] = a[1].dot(boundary_condition[z_n][0])
    
    for k in range(2, sum_):
        A[k] = np.linalg.solve(np.eye(num_entries(k, N)) - a[k].dot(A[k - 1]), b[k])
        # G[k] = np.linalg.solve(np.eye(num_entries(k, N)) - a[k].dot(A[k - 1]), (a[k].dot(G[k-1]))) 


    A[3*N - 4] = a[3*N - 4]
    # G[3*N - 4] = b[3*N - 4].dot(boundary_condition[z_n][1])
    for k in range(3*N - 5, sum_, -1):
        A[k] = np.linalg.solve(np.eye(num_entries(k, N)) - b[k].dot(A[k + 1]), a[k])
        # G[k] = np.linalg.solve(np.eye(num_entries(k, N)) - b[k].dot(A[k + 1]), (b[k].dot(G[k+1]))) 

    
    # Calculate Green's Functions
    V = {}
    k = sum_
    V[k] = np.linalg.solve(np.eye(num_entries(k, N)) - a[k].dot(A[k - 1]) - b[k].dot(A[k + 1]), 
                           C)

    for k in range(sum_ + 1, 1669):
        V[k] = A[k].dot(V[k - 1])


    toc = time.time()  
    print (toc-tic)
    return V[1668][inv_mapping[1668][(451,537,680)]]

array4 = []
# for w in np.linspace(-20, 20, 100):
array4.append (calculate(1+0.1j, N = 1000))

toc = time.time()
print ("Time taken is", toc-tic)
