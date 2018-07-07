import numpy as np
import time
import matplotlib.pyplot as plt
import time

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
        mapping[k] = {i: triples[i] for i in range(len(triples))}
        inv_mapping[k] = {triples[i]: i for i in range(len(triples))}
    return mapping, inv_mapping

def omit_sites(inv_mapping, omit):    
    site = 0
    new_inv_map = {}
    new_map = {}
    
    for sum_ in inv_mapping.keys():
        inv_sum_k_sites = {}
        sum_k_sites = {}
        site_index = 0
        for triple in inv_mapping[sum_].keys():
            if omit[site] == 1:
                inv_sum_k_sites[triple] = site_index
                sum_k_sites[site_index] = triple
                site_index += 1
            else:
                x, y, z = triple
#                 print ("I'm omitting", triple)
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
    V = {}
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
    for k in range(sum_ - 1, 0):
        V[k] = A[k].dot(V[k + 1])

    # return V[x_t + y_t + z_t]
    return V

# print(inv_mapping[(1,4,2)])
# print (coords)


# Create training data
N=20
mapping, inv_mapping = index_triple_maps(N)
X = []
Y = []
for _ in range(1000):
    x_, y_, z_ = np.random.randint(20, size=3)
    if x_ + y_ + z_ == 3*N - 3:
        continue
    omit = np.full(N**3, 1)
    omit[np.random.randint(N**3)] = 0
    _mapping, _inv_mapping, x_o, y_o, z_o = omit_sites(inv_mapping, omit)
    # Do not omit target and destination states
    if (x_ == x_o) & (y_ == y_o) & (z_ == z_o):
        continue
    

    energy = np.random.uniform(-1, 1) 
    print (x_, y_, z_, x_o, y_o, z_o, energy)
    
    true_value = (_calculate(energy+0.1j, x_, y_, z_, x_o, y_o, z_o, mapping = mapping, inv_mapping = inv_mapping))

    # Remove omitted state from true value for error calculation
    if x_o + y_o + z_o in true_value.keys():
        true_value[x_o + y_o + z_o] = np.delete(true_value[x_o + y_o + z_o], inv_mapping[x_o + y_o + z_o][(x_o, y_o, z_o)])
        true_value[x_o + y_o + z_o] = true_value[x_o + y_o + z_o][:,np.newaxis]

    
    predicted_value = (_calculate(energy+0.1j, x_, y_, z_, x_o, y_o, z_o, mapping = _mapping, inv_mapping = _inv_mapping))
    
    # Calculate error 
    for k in true_value.keys():
        error = np.absolute((predicted_value[k] - true_value[k]) / true_value[k] * 100)
        for index, triple in enumerate(_inv_mapping[k].keys()):
            (x, y, z) = triple
            x, y, z, x__, y__, z__, x_o_, y_o_, z_o_ = (np.array([x, y, z, x_, y_, z_, x_o, y_o, z_o], dtype=float) / 10.0) - 1
            X.append([x, y, z, x__, y__, z__, x_o_, y_o_, z_o_, energy])
        Y.append(error)

X = np.array(X)
print ("Error shapes", [ (i).shape for i in Y])
Y_ = np.concatenate(Y)
Y_ = (Y_ - np.mean(Y_))/np.std(Y_)
Y_.shape

from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, Y_, test_size=0.2, random_state=42)
print (X_train.shape)
from keras.models import Sequential
from keras.layers import Dense
from keras.wrappers.scikit_learn import KerasRegressor
from keras import regularizers
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import mean_squared_error


seed = 7
np.random.seed(seed)



def baseline_model():
    # create model
    model = Sequential()
    model.add(Dense(40, input_dim=10, kernel_initializer='normal', activation='relu', kernel_regularizer=regularizers.l2(0.001)))
    model.add(Dense(40, input_dim=40, kernel_initializer='normal', activation='relu', kernel_regularizer=regularizers.l2(0.001)))
    model.add(Dense(40, input_dim=40, kernel_initializer='normal', activation='relu', kernel_regularizer=regularizers.l2(0.001)))
    model.add(Dense(1, kernel_initializer='normal', kernel_regularizer=regularizers.l2(0.001)))
    # Compile model
    model.compile(loss='mean_squared_error', optimizer='adam', metrics = ["mean_squared_error"])
    return model

# estimator = KerasRegressor(build_fn=baseline_model, epochs=200, batch_size=50, verbose=1)

# # kfold = KFold(n_splits=2, random_state=seed)
# # results = cross_val_score(estimator, X_train, y_train, cv=kfold)
# # print("Results: %.2f (%.2f) MSE" % (results.mean(), results.std()))
# # print ("Results", results)

# estimator.fit(X_train, y_train)
# prediction = estimator.predict(X_test)
# print ("Mean squared error:",  mean_squared_error(y_test, prediction))
np.save("pred", X)
np.save("target", Y_)




# Reimplement since we generate a ton of data with one calculation
