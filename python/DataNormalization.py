from sklearn import preprocessing
import numpy as np

# every point is a row
# every feature is a column


# returns scaler
# to transform data need to apply scaler.transform
def DataNormalizer(unscaled_X):
    scaler = preprocessing.StandardScaler().fit(unscaled_X)
    return scaler


if __name__ =="__main__":
    X = np.array([[1., -1., 2.],
                  [2., 0., 0.],
                  [0., 1., -1.]])

    scaler = DataNormalizer(X)

    # scaled input
    X_scaled = scaler.transform(X)

    test_Point = [-1., 1., 0.]
    test_Point_scaled = scaler.transform([test_Point])

    print X_scaled