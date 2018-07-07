import numpy as np
import time
import matplotlib.pyplot as plt

e_n = 1 # energy
t_n = 1 # amplitude for hopping (|r| = 1)


N = 5 # length of cume

def gen_triples(N):
    i, j, k = (0, 0, 0)
    coords = []
    for i in range(N):
        for j in range(N):
            for k in range(N):
                coords.append((i,j,k))

    return coords


tic = time.time()   

coords = gen_triples(N)
mapping = {i: coords[i] for i in range(len(coords))}
inv_mapping = {coords[i]: i for i in range(len(coords))}
print(inv_mapping[(1,4,2)])
# print (coords)
toc = time.time()
print (toc-tic)

array = []
sites = N**3
H = np.zeros((sites, sites), dtype = complex)
for i in range(sites):
    x, y, z = mapping[i]
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


w, v = np.linalg.eig(H)
for x in np.linspace(-20, 20, 10):
    z = x + 0.1j # energy value
    
    coord_ = inv_mapping[(1, 4, 2)]
    coord = inv_mapping[(2, 2, 2)]
    G = np.sum(v[coord] * v[coord_] / (z - w))
    print (G)
    array.append (G)

print ("array", array)
plt.plot(np.linspace(-20, 20, 10), np.real(array))
plt.show()



    # for n in range(N):
    #     H[n, n] = e_n
    #     if n - 1 >= 0:
    #         H[n, n-1] = t_n
    #     if n + 1 < N:
    #         H[n, n+1] = t_n
    # w, v = np.linalg.eig(H)

    # G_ = np.zeros((N, N), dtype = complex)
    # for n in range(N):
    #     for n_ in range(N):
    #         G_[n, n_] = np.sum(v[n] * v[n_] / (z - w))

    # G__values.append(G_[100, 100])
    # print("Energy value:", z, "G: ", G_[100, 100])