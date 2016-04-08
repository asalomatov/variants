#!/mnt/xfs1/home/asalomatov/miniconda2/bin/python
import sys
sys.path.insert(0, '/mnt/xfs1/home/asalomatov/projects/variants/variants')
import train
import numpy
import pandas
import ped
import os
from sklearn.externals import joblib
from keras.models import model_from_json

m = sys.argv[1]
print sys.argv
lvl = int(sys.argv[2])
test_set_pat = sys.argv[3]
list_of_features = sys.argv[4]
ped_file = sys.argv[5]
is_keras = bool(int(sys.argv[6]))
m_pkl = joblib.load(m)

myped = ped.Ped(ped_file)
myped.addTestFile(field='ind_id', file_pat=test_set_pat)
myped.ped.dropna(subset=['test'], inplace=True)
myped.ped.reset_index(inplace=True)

test_labels = numpy.array([], dtype=int)
pred_labels = numpy.array([], dtype=int)
test_var_id = numpy.array([], dtype=str)
test_alleles = numpy.array([], dtype=str)
pred_prob = numpy.array([], dtype=float)
for i, row in myped.ped.iterrows():
    print 'processing', i, row['ind_id']
    #myped.ped.test.iat[i]
    tst = train.TrainTest(row['test'],
                          list_of_features,
                          m_pkl['y_name'],
                          m_pkl['extra_col_names'])
    if is_keras:
        tst.is_keras = True
    tst.readFeatureList()
    #tst.readTestSet()
    tst.readDataSet()
    print 'data_set shape is ', tst.data_set.shape
    tst.addLabels(level=lvl)
    print 'data_set shape is ', tst.data_set.shape
    print tst.data_set.label.value_counts()
    tst.dropNA('label')
    print 'data_set shape with non null labels is ', tst.data_set.shape
    if tst.is_keras:
        tst.model = model_from_json(m_pkl['model'])
        tst.model.load_weights(m_pkl['weights_file'])
    else:
        tst.model = m_pkl['model']
    tst.stdize = m_pkl['stdize']  # bool(int(sys.argv[6]))
    tst.trshold = m_pkl['threshold']  # float(sys.argv[7])
    tst.train_set_var_id = m_pkl['train_var_id']
    tst.data2Test()
    print 'test_set_X shape is', tst.test_set_X.shape
    tst.predictClass(tst.threshold)
    #tst.getMetrics()
    test_labels = numpy.concatenate((test_labels, tst.test_set_y))
    pred_labels = numpy.concatenate((pred_labels, tst.pred_y))
    pred_prob = numpy.concatenate((pred_prob, tst.pred_y_prob))
    test_var_id = numpy.concatenate((test_var_id, tst.test_set_var_id))
    test_alleles = numpy.concatenate((test_alleles, tst.test_set_alleles))

res = pandas.DataFrame({'test_labels': test_labels, 'pred_labels': pred_labels,
                        'pred_prob': pred_prob, 'test_var_id': test_var_id,
                        'test_var_alleles': test_alleles})
m_name = os.path.basename(m)
m_name = '.'.join(m_name.split('.')[:-1]) + '_tstlvl' + str(lvl)
res['method'] = m_name
res.to_csv(os.path.join(os.path.dirname(test_set_pat), m_name + '.csv'), index=False)

