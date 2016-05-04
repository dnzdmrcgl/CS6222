from sklearn.cross_validation import cross_val_score
from sklearn.datasets import load_iris
from sklearn.ensemble import BaggingClassifier
from sklearn.neighbors import KNeighborsClassifier

bagging = BaggingClassifier(KNeighborsClassifier(), max_samples=0.5, max_features=0.5)
iris = load_iris()

bagging.fit(iris.data, iris.target)
print iris.target

scores = cross_val_score(bagging, iris.data, iris.target)
print scores.mean()