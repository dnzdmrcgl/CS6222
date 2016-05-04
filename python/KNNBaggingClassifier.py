from sklearn.cross_validation import cross_val_score
from sklearn.datasets import load_iris
from sklearn.ensemble import BaggingClassifier
from sklearn.neighbors import KNeighborsClassifier



def GetAdaKNNBaggingClassifier(X, y):
    clf  = BaggingClassifier(KNeighborsClassifier(), max_samples=0.5, max_features=0.5)
    clf.fit(X, y)
    return clf

def EvaluateKNNBagging(X, y, scoring = None):
    clf  = BaggingClassifier(KNeighborsClassifier(), max_samples=0.5, max_features=0.5)
    scores = cross_val_score(clf, X, y, scoring= scoring)
    return scores.mean()

if __name__ == "__main__":
    iris = load_iris()
    print EvaluateKNNBagging(iris.data, iris.target)