#!/mnt/xfs1/home/asalomatov/miniconda2/bin/python
from __future__ import print_function
import sys
#sys.path.insert(0, '/mnt/xfs1/home/asalomatov/projects/update_vars/variants/variants')
import train
import numpy
import pandas
import os
#from ggplot import *
from sklearn.externals import joblib
from keras.models import model_from_json

print(sys.argv)

m = sys.argv[1]
lvl = int(sys.argv[2])
n_extra = int(sys.argv[3])
prob_threshold = float(sys.argv[4])
known_vars = sys.argv[5]
extra_vars = sys.argv[6]
is_keras = bool(int(sys.argv[7]))
m_pkl = joblib.load(m)
list_of_features = m_pkl['features']


tst = train.TrainTest(known_vars,
                      list_of_features,
                      m_pkl['y_name'],
                      m_pkl['extra_col_names'])
# +
#                      ['DP_offspring', 'DP_father', 'DP_mother'])
tst.feature_list = list_of_features
if is_keras:
    tst.is_keras = True
tst.readDataSet()
tst.addLabels(level=lvl)

print('data_set shape is %s' % ' '.join(map(str, tst.data_set.shape)))
if tst.data_set.empty:
    sys.exit('data set is empty')

#n_extra = tst.data_set.shape[0] #roughly balanced classes
print('adding %s extra negative examples' % n_extra)
if n_extra > 0:
    tst.readExtraVars(extra_vars, n_extra=n_extra)
tst.dropNA('label')
print('data_set shape is %s' % ' '.join(map(str, tst.data_set.shape)))
print('label balance is ')
print(tst.data_set.label.value_counts())
if tst.is_keras:
    tst.model = model_from_json(m_pkl['model'])
    tst.model.load_weights(m_pkl['weights_file'])
else:
    tst.model = m_pkl['model']
tst.stdize = m_pkl['stdize']  # bool(int(sys.argv[6]))
tst.threshold = prob_threshold
print('probability treshold is %s' % tst.threshold)
tst.train_set_var_id = [] #m_pkl['train_var_id']
tst.data2Test()
print('test_set_X shape is')
print(tst.test_set_X.shape)
tst.predictClass(tst.threshold)
tst.getMetrics()
tst.perf_mertics['method'] = tst.method
tst.perf_mertics['prob_cutoff'] = tst.threshold


# myplot = ggplot(tst.roc, aes(x='fpr', y='tpr')) +\
#     geom_line() +\
#     geom_abline(linetype='dashed')
# myplot1 = ggplot(tst.roc, aes('threshold')) +\
#     geom_line(aes(y='tpr')) +\
#     geom_line(aes(y='fpr'))
# ggsave(plot=myplot, filename='roc_curve_1.png')
# ggsave(plot=myplot1, filename='roc_curve_2.png')


#m_name = os.path.basename(m)
#m_name = '.'.join(m_name.split('.')[:-1]) + '_tstlvl' + str(lvl)
#res['method'] = m_name
#res.to_csv(os.path.join(os.path.dirname(test_set_pat), m_name + '.csv'), index=False)

