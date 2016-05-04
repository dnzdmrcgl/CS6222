import cStringIO
import numpy as np
from sklearn.cross_validation import train_test_split
X, y = np.arange(10).reshape((5, 2)), range(5)


# Splits dataset into test and train parts
def SplitTestTrain(X_total, y_total):
    X_train, X_test, y_train, y_test = train_test_split(X_total, y_total, test_size=0.25, random_state=42)
    return X_train, y_train, X_test, y_test



#print  float(False)
with open('../features/real1Features.csv', "r") as myfile:
    data = myfile.read().replace('TRUE', '1.0').replace('FALSE', '0.0')
input_data = np.genfromtxt(cStringIO.StringIO(data), delimiter=",", skip_header=1, case_sensitive=False)
X_full = input_data[:, :-1]
#print X_full.shape
y_full = input_data[:, -1]
#print y_full.shape
