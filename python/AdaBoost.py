from sklearn.cross_validation import cross_val_score
from sklearn.datasets import load_iris
from sklearn.ensemble import AdaBoostClassifier


# AdaBoost classifier for ensemble learning
# doesn't need normalization


def GetAdaBoostClassifier(X, y):
    clf = AdaBoostClassifier(n_estimators=100)
    clf.fit(X, y)
    return clf

# Evaluate by cross validation
def EvaluateAdaBoost(X, y, scoring = None):
    clf = AdaBoostClassifier(n_estimators=1000)
    scores = cross_val_score(clf, X, y, scoring= scoring)
    return scores.mean()

if __name__ == "__main__":
    iris = load_iris()
    clf = AdaBoostClassifier(n_estimators=100)
    scores = cross_val_score(clf, iris.data, iris.target)
    print scores.mean()
