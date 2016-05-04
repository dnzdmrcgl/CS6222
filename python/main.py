from sklearn.datasets import load_iris

from AdaBoost import EvaluateAdaBoost
from RandomForestWithFeatureExtraction import EvaluateRF
from KNNBaggingClassifier import EvaluateKNNBagging
from DataNormalization import DataNormalizer
from utils import DataSetReader

def testWithFileName(filename):
    print filename
    X, y = DataSetReader(filename)

    # scoring = None
    scoring = 'f1'

    # iris = load_iris()
    # X, y = iris.data, iris.target

    print "Ada Boost score " + str(EvaluateAdaBoost(X, y, scoring))

    # Other two classifiers depend on distance, hence need to normalize data
    scaler = DataNormalizer(X)
    scaled_X = scaler.transform(X)
    print "Random forest with feature extraction score " + str(EvaluateRF(scaled_X, y, scoring))

    print "KNN Bagging score " + str(EvaluateKNNBagging(scaled_X, y, scoring))
    print

if __name__ == '__main__':
    filename = '../features/real1Features.csv'
    testWithFileName(filename)
    for i in range(1,6):
        filename = '../features/syn' + str(i) + 'Features.csv'
        testWithFileName(filename)