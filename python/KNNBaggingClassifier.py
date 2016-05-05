from sklearn.cross_validation import cross_val_score
from sklearn.datasets import load_iris
from sklearn.ensemble import BaggingClassifier
from sklearn.neighbors import KNeighborsClassifier

# Bagging method for ensemble learning
# As a weak classifier we use kNN - k nearest neighbours


def GetAdaKNNBaggingClassifier(X, y):
    clf  = BaggingClassifier(KNeighborsClassifier(), max_samples=0.7, max_features=0.7)
    clf.fit(X, y)
    return clf

# Evaluate by cross validation
def EvaluateKNNBagging(X, y, scoring = None):
    clf  = BaggingClassifier(KNeighborsClassifier(), max_samples=0.7, max_features=0.7)
    scores = cross_val_score(clf, X, y, scoring= scoring)
    return scores.mean()

if __name__ == "__main__":
    iris = load_iris()
    print EvaluateKNNBagging(iris.data, iris.target)