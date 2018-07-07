import numpy as np
X = np.load("pred.npy")
Y_ = np.load("target.npy")


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
    model.add(Dense(100, input_dim=10, kernel_initializer='normal', activation='relu', kernel_regularizer=regularizers.l2(0.0001)))
    model.add(Dense(100, input_dim=100, kernel_initializer='normal', activation='relu', kernel_regularizer=regularizers.l2(0.0001)))
    model.add(Dense(100, input_dim=100, kernel_initializer='normal', activation='relu', kernel_regularizer=regularizers.l2(0.0001)))
    model.add(Dense(1, kernel_initializer='normal', kernel_regularizer=regularizers.l2(0.0001)))
    # Compile model
    model.compile(loss='mean_squared_error', optimizer='adam', metrics = ["mean_squared_error"])
    return model

estimator = KerasRegressor(build_fn=baseline_model, epochs=2000, batch_size=50, verbose=1)

# kfold = KFold(n_splits=2, random_state=seed)
# results = cross_val_score(estimator, X_train, y_train, cv=kfold)
# print("Results: %.2f (%.2f) MSE" % (results.mean(), results.std()))
# print ("Results", results)

estimator.fit(X_train, y_train)
prediction = estimator.predict(X_test)
print ("Mean squared error:",  mean_squared_error(y_test, prediction))
