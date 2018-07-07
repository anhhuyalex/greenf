import numpy as np
import time
import matplotlib.pyplot as plt

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






def calculate6_8_9_g_3_4_7(z_n, t_n = 1.0, e_n = 1.0, N = 20, k = 1, c = complex):

    mapping, inv_mapping = index_triple_maps(N)

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
                if x + 1 <= N - 1:
                    b_k[index, inv_mapping[k + 1][x + 1, y, z]] = amplitude

                if y + 1 <= N - 1:
                    b_k[index, inv_mapping[k + 1][x, y + 1, z]] = amplitude

                if z + 1 <= N - 1:
                    b_k[index, inv_mapping[k + 1][x, y, z + 1]] = amplitude

            # Calculate a_k
            if 0 < k:
                if x - 1 >= 0:
                    a_k[index, inv_mapping[k - 1][x - 1, y, z]] = amplitude

                if y - 1 >= 0:
                    a_k[index, inv_mapping[k - 1][x, y - 1, z]] = amplitude

                if z - 1 >= 0:
                    a_k[index, inv_mapping[k - 1][x, y, z - 1]] = amplitude

        a[k] = a_k
        b[k] = b_k

    # Calculate C
    x_, y_, z_ = (1, 4, 2)
    sum_ = sum([x_, y_, z_])
    C = np.zeros(shape = num_entries(sum_, N), dtype = complex)
    C[inv_mapping[sum_][(x_, y_, z_)]] = 1. / (z_n - e_n)
    # print (C)


    # Calculate Multiplicative Factors 
    A = {}
    A[0] = b[0]
    for k in range(1, sum_):
        A[k] = np.linalg.solve(np.eye(num_entries(k, N)) - a[k].dot(A[k - 1]), b[k])
        # print (A[k])


    A[3*N - 3] = a[3*N - 3]
    # print (a[3*N - 3])
    for k in range(3*N - 4, sum_, -1):
        A[k] = np.linalg.solve(np.eye(num_entries(k, N)) - b[k].dot(A[k + 1]), a[k])

    # Calculate Green's Functions
    V = {}
    k = sum_
    V[k] = np.linalg.solve(np.eye(num_entries(k, N)) - a[k].dot(A[k - 1]) - b[k].dot(A[k + 1]), C)
    # print(k, V[k])
    # for k in range(sum_ + 1, 3*N - 2):
    for k in range(sum_ - 1, 5, -1):
        V[k] = A[k].dot(V[k + 1])
        # print(k, V[k])


    toc = time.time()  
    print (toc-tic)
    return V[6][inv_mapping[6][(2, 2, 2)]]

array = []
for w in np.linspace(-20, 20, 10):
    array.append (calculate6_8_9_g_3_4_7(w+0.1j, N = 5))

plt.plot(np.linspace(-20, 20, 10), np.real(array))
plt.show()
# print (inv_mapping[sum_])
# print (inv_mapping[sum_][(x_, y_, z_)])



# N = range(2, 200, 10)
# my_time = [] 
# for n in N:
#     tic = time.time()
#     for k in range(3*n - 2):
#         generate_triples(n, k)
#     toc = time.time()
#     print (n, toc-tic)
#     my_time.append(toc-tic)

# plt.plot(N, my_time)
# plt.show()
# print (generate_triples(10, 18))

