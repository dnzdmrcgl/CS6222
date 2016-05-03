# NOTE: should do it before normalization

# feature is a column,
# stolbetz
import numpy as np
from sklearn.preprocessing import Imputer
imp = Imputer(missing_values='NaN', strategy='mean', axis=0)
X_full = np.array([[1, 2], [np.nan, 3], [7, 6]])
imp.fit(X_full)

X_tests= [[np.nan, 2], [6, np.nan]]

print(imp.transform(X_tests))
