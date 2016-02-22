import sys
sys.path.insert(0, '/mnt/xfs1/home/asalomatov/projects/variants/variants')
import train
import numpy
import pandas
import ped
from sklearn.externals import joblib

model_list = ['GBM_depth1_lvl1_10000.pkl', 'GBM_depth1_lvl2_10000.pkl', 'GBM_depth1_lvl3_10000.pkl', 'GBM_depth1_lvl4_10000.pkl', 'LinearSVM_C1_lvl1_10000.pkl', 'LinearSVM_C1_lvl2_10000.pkl', 'LinearSVM_C1_lvl3_10000.pkl', 'LinearSVM_C1_lvl4_10000.pkl', 'LogReg_10000.pkl', 'polySVM_deg3_C1_lvl1_10000.pkl', 'polySVM_deg3_C1_lvl2_10000.pkl', 'polySVM_deg3_C1_lvl3_10000.pkl', 'polySVM_deg3_C1_lvl4_10000.pkl', 'RandForest_depth1_lvl1_10000.pkl', 'RandForest_depth1_lvl2_10000.pkl', 'RandForest_depth1_lvl3_10000.pkl', 'RandForest_depth1_lvl4_10000.pkl', 'rbfSVM_g07_C1_lvl1_10000.pkl', 'rbfSVM_g07_C1_lvl2_10000.pkl', 'rbfSVM_g07_C1_lvl3_10000.pkl', 'rbfSVM_g07_C1_lvl4_10000.pkl']

jobB
test_labels = numpy.array([], dtype=int)
pred_labels = numpy.array([], dtype=int)
test_var_id = numpy.array([], dtype=str)
pred_prob = numpy.array([], dtype=float)
for i, row in myped.ped.iterrows():
    print 'processing', i, row['ind_id']
    #myped.ped.test.iat[i]

    tst = train.TrainTest(row['test'],
                    '/mnt/xfs1/home/asalomatov/projects/variants/variants/ssc_wes_features.txt',
                    'status')
    tst.readFeatureList()
    tst.readDataSet()
    print 'data_set shape is ', tst.data_set.shape
    tst.addLabels(level=4)
    print 'data_set shape is ', tst.data_set.shape
    tst.data_set.label.isnull().sum()
    print tst.data_set.label.value_counts()
    #estGB_classifier = joblib.load('estGB_classifier.pkl')
    tst.model = trn.model
    tst.train_set_var_id = trn.train_set_var_id
    tst.data2Test()
    print 'test_set_X shape is', tst.test_set_X.shape
    tst.predictClass(threshold=0.5)
    #tst.getMetrics()
    test_labels = numpy.concatenate((test_labels, tst.test_set_y))
    pred_labels = numpy.concatenate((pred_labels, tst.pred_y))
    pred_prob = numpy.concatenate((pred_prob, tst.pred_y_prob))
    test_var_id = numpy.concatenate((test_var_id, tst.test_set_var_id))

