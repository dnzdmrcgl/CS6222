
# A class to treat missing values  in the data - they are replaced by the mean of the corresponding feature


# feature is a column,
# stolbetz
import numpy as np
from sklearn.preprocessing import Imputer

# Return Imputer instance
# to fill missing values need to apply transform
def MisingValuesFiller(X_train):
    imp = Imputer(missing_values='NaN', strategy='mean', axis=0)
    X_full = np.array(X_train)
    imp.fit(X_full)
    #print X_full
    X_tests= [[np.nan, 2], [6, np.nan]]

    #print(imp.transform(X_tests))
    return imp

if __name__ == "__main__":
    X_train = [[1, 2], [np.nan, 3], [7, 6]]

    imp = MisingValuesFiller(X_train)
    print imp.transform(X_train)