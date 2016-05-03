from sklearn import preprocessing
import numpy as np

# every point is a row
# every feature is a column
X = np.array([[ 1., -1.,  2.],
               [ 2.,  0.,  0.],
               [ 0.,  1., -1.]])

scaler = preprocessing.StandardScaler().fit(X)

# scaled input
X_scaled = scaler.transform(X)

test_Point = [-1., 1., 0.]
test_Point_scaled = scaler.transform([test_Point])

print X_scaled


