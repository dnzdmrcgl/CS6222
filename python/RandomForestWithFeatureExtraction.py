from sklearn.cross_validation import cross_val_score
from sklearn.datasets import load_iris
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.pipeline import Pipeline
from sklearn.svm import LinearSVC



def GetRFClassifier(X, y):
  clf = Pipeline([
    ('feature_selection', SelectFromModel(LinearSVC(penalty="l1", dual=False))),
    ('classification', RandomForestClassifier())
  ])
  clf.fit(X, y)
  return clf


def EvaluateRF(X, y, scoring = None):
    clf = Pipeline([
      ('feature_selection', SelectFromModel(LinearSVC(penalty="l1", dual=False))),
      ('classification', RandomForestClassifier())
    ])
    scores = cross_val_score(clf, X, y, scoring= scoring)
    return scores.mean()

if __name__ == '__main__':
  iris = load_iris()
  X, y = iris.data, iris.target
  GetRFClassifier(X, y)
  print EvaluateRF(X, y)
