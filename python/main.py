from sklearn.datasets import load_iris

from AdaBoost import EvaluateAdaBoost
from RandomForestWithFeatureExtraction import EvaluateRF
from KNNBaggingClassifier import EvaluateKNNBagging
from DataNormalization import DataNormalizer
from utils import DataSetReader

def testWithFileName(input_data_filename,output_file, scoring):
    print input_data_filename
    output_path.write(input_data_filename + '\n')
    X, y = DataSetReader(input_data_filename)

    # scoring = None
    #scoring = 'f1'

    # iris = load_iris()
    # X, y = iris.data, iris.target

    s1 = EvaluateAdaBoost(X, y, scoring)
    print "Ada Boost score " + str(s1)
    output_file.write("Ada Boost score " + str(s1) + '\n')
    # Other two classifiers depend on distance, hence need to normalize data
    scaler = DataNormalizer(X)
    scaled_X = scaler.transform(X)

    s2 = EvaluateRF(scaled_X, y, scoring)
    print "Random forest with feature extraction score " + str(s2)
    output_file.write( "Random forest with feature extraction score " + str(s2) + '\n')

    s3 = EvaluateKNNBagging(scaled_X, y, scoring)
    print "KNN Bagging score " + str(s3)
    output_file.write( "KNN Bagging score " + str(s3) + '\n\n')
    print

if __name__ == '__main__':
    output_path = 'predcitions_scores.txt'
    myfile = open(output_path,'w')
    scores = [ 'recall', 'precision', 'f1']
    #scores = [ 'recall']
    for score in scores:
        print "score is " + score
        myfile.write("score is " + score + '\n')


        for i in range(1,3):
            filename = '../features/real' + str(i) +'Features.csv'
            testWithFileName(filename,myfile, score)

        for i in range(1,6):
            filename = '../features/syn' + str(i) + 'Features.csv'
            testWithFileName(filename,myfile, score)

        print
        print
    myfile.close()
    print "Done"