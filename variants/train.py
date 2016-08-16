from __future__ import print_function
import sys
#sys.path.insert(0, '/mnt/xfs1/home/asalomatov/projects/variants/variants')
import ped
import variants
import func
import pandas
import numpy
import os
import features
import train
from multiprocessing import Pool
from sklearn.cross_validation import train_test_split
from sklearn.ensemble import (RandomForestClassifier,
                              GradientBoostingClassifier,
                              GradientBoostingRegressor)
from sklearn.linear_model import LogisticRegression
from sklearn import svm, metrics
from sklearn.externals import joblib
from sklearn.preprocessing import binarize
from sklearn.preprocessing import scale
from sklearn.decomposition import PCA
#from unbalanced_dataset import SMOTE
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#import seaborn as sns
from keras.models import Sequential, model_from_json
from keras.layers.core import Dense, Activation, Dropout
import yaml
import argparse

# %matplotlib qt
# sns.set()
# plt.rcParams['font.size'] = 14
# almost_black = '#262626'


class TrainTest:

    def __init__(self, data_set_file, feature_list_file, y_name=[],
                 extra_col_names=[]):
        self.data_set = None
        self.train_set_X = None
        self.train_set_y = None
        self.train_set_var_id = []
        self.test_set_var_id = []
        self.train_set_alleles = []
        self.test_set_alleles = []
        self.test_set_DP_father = []
        self.test_set_DP_mother = []
        self.test_set_X = None
        self.test_set_y = None
        self.pred_y = None
        self.pred_y_prob = None
        self.feature_list = None
        self.feature_list_file = feature_list_file
        self.data_set_file = data_set_file
        self.y_name = y_name # a list of length 0 or 1
        self.extra_column_names = extra_col_names  # a list of length >= 0
        self.method = None
        self.model = None
        self.perf_mertics = None
        self.roc = None
        self.stdize = True
        self.feature_importance = None
        self.threshold = 0.5
        self.is_keras = False

    def addAlleleBalance(self, mydf=None):
        if mydf is None:
            self.data_set['allele_balance_offspring'] = self.data_set['REF_count_offspring']/self.data_set['DP_offspring']
            self.data_set['allele_balance_father'] = self.data_set['REF_count_father']/self.data_set['DP_father']
            self.data_set['allele_balance_mother'] = self.data_set['REF_count_mother']/self.data_set['DP_mother']
        else:
            mydf['allele_balance_offspring'] = mydf['REF_count_offspring']/mydf['DP_offspring']
            mydf['allele_balance_father'] = mydf['REF_count_father']/mydf['DP_father']
            mydf['allele_balance_mother'] = mydf['REF_count_mother']/mydf['DP_mother']
            return mydf

    def addLabels(self, level):
        """level > 1 are specific to SF data.
        level == 1 applies to any set. In column self.y_name positive 
        examples should be denoted with 'Y', and negative ones with 'N'
        level = 1 - Y vs N
        level = 2 - Y + ND & both vs N
        level = 3 - Y + ND vs N
        level = 4 - Y + ND vs all other
        """
        if level == 0:
            self.data_set['label'] = -1
        elif level == 1:
            self.data_set['label'] = None
            self.data_set.ix[self.data_set[self.y_name[0]].isin(['Y']), 'label'] = 1
            self.data_set.ix[self.data_set[self.y_name[0]].isin(['N']), 'label'] = 0
        elif level == 2:
            self.data_set['label'] = None
            self.data_set.ix[self.data_set.status.isin(['Y']), 'label'] = 1.
            c1 = self.data_set.status.isin(['ND'])
            c2 = self.data_set.descr.isin(['both'])
            self.data_set.ix[c1 & c2, 'label'] = 1
            self.data_set.ix[self.data_set.status.isin(['N']), 'label'] = 0
        elif level == 3:
            self.data_set['label'] = None
            self.data_set.ix[self.data_set.status.isin(['Y']), 'label'] = 1.
            c1 = self.data_set.status.isin(['ND'])
            self.data_set.ix[c1, 'label'] = 1
            self.data_set.ix[self.data_set.status.isin(['N']), 'label'] = 0
        elif level == 4:
            self.data_set['label'] = 0
            c1 = self.data_set.status.isin(['ND', 'Y'])
            self.data_set.ix[c1, 'label'] = 1
        elif level == 5:
            self.data_set['label'] = 0
            self.data_set.ix[self.data_set.status.isin(['Y']), 'label'] = 1.
            c1 = self.data_set.status.isin(['ND'])
            c2 = self.data_set.descr.isin(['both'])
            self.data_set.ix[c1 & c2, 'label'] = 1
        elif level == 6:
            self.data_set['label'] = 0
            self.data_set.ix[self.data_set.status.isin(['Y']), 'label'] = 1.
            c1 = self.data_set.status.isin(['ND'])
            c2 = self.data_set.descr.isin(['both'])
            self.data_set.ix[c1 & c2, 'label'] = 1
            self.data_set.ix[c1 & (~c2), 'label'] = None
        elif level == 7:
            self.data_set['label'] = 0
            self.data_set.ix[self.data_set.status.isin(['Y']), 'label'] = 1.
            c1 = self.data_set.status.isin(['ND'])
            c2 = self.data_set.descr.isin(['both'])
            self.data_set.ix[c1 & c2, 'label'] = 1
        elif level == 8:
            self.data_set['label'] = 0
            c_Y = self.data_set[self.y_name[0]].isin(['Y'])
            c_N = self.data_set[self.y_name[0]].isin(['N'])
            self.data_set.ix[c_Y, 'label'] = 1
            self.data_set.ix[c_N, 'label'] = 0
            c_ND = self.data_set.status.isin(['ND'])
            self.data_set.ix[c_ND, 'label'] = None
        elif level == 9:
            self.data_set['label'] = 0
            c_Y = self.data_set[self.y_name[0]].isin(['Y'])
            c_N = self.data_set[self.y_name[0]].isin(['N'])
            c_Krumm = self.data_set.descr.isin(['Krumm'])
            c_both = self.data_set.descr.isin(['both'])
            c_ioss = self.data_set.descr.isin(['Iossifov'])
            self.data_set.ix[c_Y & (c_Krumm | c_both), 'label'] = 1
            self.data_set.ix[c_N & (c_Krumm | c_both), 'label'] = 0
            c_ND = self.data_set.status.isin(['ND'])
            self.data_set.ix[c_ND | c_ioss, 'label'] = None
        elif level == 10:  # columbia data only 
            self.data_set['label'] = 0
            c_Y = self.data_set[self.y_name[0]].isin(['Y'])
            c_N = self.data_set[self.y_name[0]].isin(['N'])
            c_after = self.data_set.descr.isin(['after'])
            c_before = self.data_set.descr.isin(['before'])
            self.data_set.ix[c_Y & c_after, 'label'] = 1
            self.data_set.ix[c_N & c_after, 'label'] = 0
            self.data_set.ix[c_before, 'label'] = None
        else:
            sys.exit('Unknown level, exiting...\n')

    def dropNA(self, field):
        self.data_set = self.data_set[~self.data_set[field].isnull()]

    def addVarID(self, df):
        df['var_id'] = df.ind_id.astype(str) +\
                                 '_' +\
                                 df.CHROM.astype(str) +\
                                 '_' +\
                                 df.POS.astype(str)

    def addAlleles(self, df):
        df['offspring_alleles'] = df.REF_base_offspring +\
                                 '_' +\
                                 df.REF_count_offspring.astype(int).astype(str) +\
                                 '_' +\
                                 df.ALT_base_offspring.astype(str)

    def addAllelesBalByDP(self, df):
        df['allele_balance_by_DP_offspring'] = df['allele_balance_offspring'] /\
                                               df['DP_offspring'].astype(float)
        df['allele_balance_by_DP_father'] = df['allele_balance_father'] /\
                df['DP_father'].astype(float)
        df['allele_balance_by_DP_mother'] = df['allele_balance_mother'] /\
                            df['DP_mother'].astype(float)

    def readFeatureList(self):
        """features_file contains names of the columns to be
        used as features, one per line
        """
        with open(self.feature_list_file, 'r') as ff:
            self.feature_list = ff.readlines()
            self.feature_list = [x.strip('\n') for x in self.feature_list]

    def readDataSet(self):
        self.data_set = pandas.read_csv(self.data_set_file, sep='\t')
        self.addAlleleBalance()
        self.addAllelesBalByDP(self.data_set)
        self.addVarID(self.data_set)
        self.addAlleles(self.data_set)
        self.data_set = self.data_set[self.feature_list +
                                      self.y_name +
                                      ['var_id'] + ['offspring_alleles'] +
                                      self.extra_column_names]

    def addExtraColumns(self):
        self.addAlleleBalance()
        self.addAllelesBalByDP(self.data_set)
        self.addVarID(self.data_set)
        self.addAlleles(self.data_set)
        self.data_set = self.data_set[self.feature_list +
                                      self.y_name +
                                      ['var_id'] + ['offspring_alleles'] +
                                      self.extra_column_names]

    def readTestSet(self):
        self.data_set = pandas.read_csv(self.data_set_file, sep='\t')
        self.addAlleleBalance()
        self.addAllelesBalByDP(self.data_set)
        self.addVarID(self.data_set)
        self.addAlleles(self.data_set)
        self.data_set = self.data_set[self.feature_list + ['var_id'] +
                                      ['offspring_alleles']]

    def readExtraVars(self, file_name,  n_extra=0):
        if n_extra > 0:
            x = pandas.read_csv(file_name, sep='\t', nrows=n_extra)
            x = x[x.status.isnull()]
            x = self.addAlleleBalance(x)
            self.addAllelesBalByDP(x)
            self.addVarID(x)
            self.addAlleles(x)
            x = x[self.feature_list + self.y_name + ['var_id'] +
                  ['offspring_alleles'] + self.extra_column_names]
            self.data_set = pandas.concat([self.data_set, x], axis=0)

    def splitTrainTest(self, trn_size=0.5, 
                       rnd_state=0, over_sample=''):
        """over_sample one of [None, 'SMOT', 'SMOT_bl1', 'SMOT_bl2', 'SMOT_svm']
        """
        X = self.data_set[self.feature_list + ['offspring_alleles'] + ['var_id']]
        y = self.data_set.label.astype(int)
        X_tr, X_te, y_tr, y_te = train_test_split(X, y, train_size=trn_size,
                                                  random_state=rnd_state,
                                                  stratify=y)
        self.train_set_X = X_tr[self.feature_list].values
        self.test_set_X = X_te[self.feature_list].values
        self.train_set_y = y_tr.values
        self.test_set_y = y_te.values
        self.train_set_var_id = list(X_tr.var_id)
        self.test_set_var_id = list(X_te.var_id)
        self.train_set_alleles = list(X_tr.offspring_alleles)
        self.test_set_alleles = list(X_te.offspring_alleles)
        if self.stdize:
            self.train_set_X = scale(self.train_set_X)
            self.test_set_X = scale(self.test_set_X)
#        if over_sample is not '':
#            ratio = float(numpy.count_nonzero(self.train_set_y == 1)) /\
#                    float(numpy.count_nonzero(self.train_set_y == 0))
#            smote = None
#            verbose = True
#            if ratio < 1:
#                ratio = 1 / ratio
#            print('upsample ratio: %s' % ratio)
#            print('smote type: %s' % over_sample)
#            if over_sample == 'SMOT':
#                smote = SMOTE(ratio=ratio, verbose=verbose, kind='regular')
#            elif over_sample == 'SMOT_bl1':
#                smote = SMOTE(ratio=ratio, verbose=verbose, kind='boderline1')
#            elif over_sample == 'SMOT_bl2':
#                smote = SMOTE(ratio=ratio, verbose=verbose, kind='boderline2')
#            elif over_sample == 'SMOT_svm':
#                smote = SMOTE(ratio=ratio, verbose=verbose, kind='svm')
#            else:
#                print('unknown SMOTE type')
#            if smote is not None:
#                print('applying SMOTE...')
#                self.train_set_X, self.train_set_y = smote.fit_transform(self.train_set_X,
#                                                                         self.train_set_y)

    def data2Test(self):
        c1 = self.data_set.var_id.isin(self.train_set_var_id)
        print(str(sum(c1)) + ' variants will be removed')
        self.test_set_X = self.data_set[self.feature_list][~c1].values
        self.test_set_y = self.data_set.label[~c1].astype(int).values
        self.test_set_var_id = self.data_set.var_id[~c1].astype(str).values
        self.test_set_alleles = self.data_set.offspring_alleles[~c1].astype(str).values
        self.test_set_DP_offspring = self.data_set.DP_offspring[~c1].astype(int).values
        self.test_set_DP_father = self.data_set.DP_father[~c1].astype(int).values
        self.test_set_DP_mother = self.data_set.DP_mother[~c1].astype(int).values
        if self.stdize:
            print('scaling data')
            self.test_set_X = scale(self.test_set_X)

    def keepPosOnly(self):
        c1 = self.train_set_y == 1
        self.train_set_alleles = self.train_set_alleles[c1]
        self.train_set_var_id = self.train_set_var_id[c1]
        self.train_set_X = self.train_set_X[c1]
        self.train_set_y = self.train_set_y[c1]

    def fitClassifier(self):
        # fit estimator
        self.model.fit(self.train_set_X, self.train_set_y)
        if hasattr(self.model, 'feature_importances_'):
            self.feature_importance = pandas.DataFrame({
                'contrib': self.model.feature_importances_,
                'name': self.feature_list})
            self.feature_importance.sort_values(['contrib'],
                                                ascending=[False],
                                                inplace=True)
            self.feature_importance.to_csv(self.method +
                                           '_feature_contrib.csv',
                                           index=False)
        #model_descr = {'model': self.model,
        #               'train_var_id': self.train_set_var_id,
        #               'stdize': self.stdize,
        #               'features': self.feature_list,
        #               'feature_importance': self.feature_importance,
        #               'y_name': self.y_name,
        #               'extra_col_names': self.extra_column_names,
        #               'method': self.method,
        #               'threshold': self.threshold}
        #joblib.dump(model_descr, self.method + '.pkl')

    def fitClassifierOneClass(self):
        # fit estimator
        # train on positive labels only
        self.keepPosOnly()
        self.model.fit(self.train_set_X, self.train_set_y)
        if hasattr(self.model, 'feature_importances_'):
            self.feature_importance = pandas.DataFrame({
                'contrib': self.model.feature_importances_,
                'name': self.feature_list})
            self.feature_importance.sort(['contrib'], ascending=[False], inplace=True)
            self.feature_importance.to_csv(self.method + '_feature_contrib.csv',
                                           index=False,
                                           sep='\t')

    def pklModel(self, output_dir='./'):
        func.runInShell('mkdir -p ' + output_dir)
        model_descr = {'model': self.model,
                       'train_var_id': self.train_set_var_id,
                       'stdize': self.stdize,
                       'features': self.feature_list,
                       'feature_importance': self.feature_importance,
                       'y_name': self.y_name,
                       'extra_col_names': self.extra_column_names,
                       'method': self.method,
                       'threshold': self.threshold,
                       'metrics': self.perf_mertics,
                       'roc': self.roc,
                       'is_keras': 0}
        joblib.dump(model_descr, os.path.join(output_dir, self.method + '.pkl'))

    def saveKerasModel(self, output_dir='./'):
        weights_file = os.path.join(output_dir, self.method + '.h5')
        self.model.save_weights(weights_file, overwrite=True)
        json_data = self.model.to_json()
        model_descr = {'model': json_data,
                       'weights_file': weights_file,
                       'train_var_id': self.train_set_var_id,
                       'stdize': self.stdize,
                       'features': self.feature_list,
                       'feature_importance': self.feature_importance,
                       'y_name': self.y_name,
                       'extra_col_names': self.extra_column_names,
                       'method': self.method,
                       'threshold': self.threshold,
                       'metrics': self.perf_mertics,
                       'is_keras': 1}
        joblib.dump(model_descr, os.path.join(output_dir, self.method + '.pkl'))

    def predictClass(self, threshold=0.5):
        # prediction
        #self.pred_y = self.model.predict(self.test_set_X)
        if self.is_keras:
            self.pred_y_prob = self.model.predict_proba(self.test_set_X)[:, 0]
        else:
            self.pred_y_prob = self.model.predict_proba(self.test_set_X)[:, 1]
        self.pred_y = binarize(self.pred_y_prob.reshape(1, -1),
                               threshold)[0].astype(int)

    def predictClassOneClass(self):
        # prediction
        self.pred_y = self.model.predict(self.test_set_X).astype(int)
        self.pred_y[self.pred_y == -1] = 0

    def getMetrics(self):
        """output various classification metrics"""
        confusion = metrics.confusion_matrix(self.test_set_y, self.pred_y)
        TP = confusion[1, 1]
        TN = confusion[0, 0]
        FP = confusion[0, 1]
        FN = confusion[1, 0]
        accuracy = metrics.accuracy_score(self.test_set_y, self.pred_y)
        null_accuracy = max(self.test_set_y.mean(), 1 - self.test_set_y.mean())
        sensitivity = metrics.recall_score(self.test_set_y, self.pred_y)
        specificity = TN / float(FP + TN)
        false_pos_rate = FP / float(FP + TN)
        precision = TP / float(TP + FP)
        f1_score = metrics.f1_score(self.test_set_y, self.pred_y)
        self.perf_mertics = pandas.DataFrame()
        self.perf_mertics['TP'] = [TP]
        self.perf_mertics['TN'] = [TN]
        self.perf_mertics['FP'] = [FP]
        self.perf_mertics['FN'] = [FN]
        self.perf_mertics['null_accuracy'] = [null_accuracy]
        self.perf_mertics['accuracy'] = [accuracy]
        self.perf_mertics['sensitivity'] = [sensitivity]
        self.perf_mertics['specificity'] = [specificity]
        self.perf_mertics['false_pos_rate'] = [false_pos_rate]
        self.perf_mertics['precision'] = [precision]
        self.perf_mertics['f1_score'] = [f1_score]
        if len(numpy.unique(self.test_set_y)) > 1 and self.pred_y_prob is not None:
            roc_auc = metrics.roc_auc_score(self.test_set_y, self.pred_y_prob)
            self.perf_mertics['roc_auc'] = [roc_auc]
            fpr, tpr, thresholds = metrics.roc_curve(self.test_set_y,
                                                     self.pred_y_prob)
            roc_df = pandas.DataFrame({'fpr': fpr,
                                       'tpr': tpr,
                                       'threshold': thresholds})
            self.roc = roc_df
        print(self.perf_mertics)
        #print metrics.classification_report(y_test, y_pred)
        #plt.plot(fpr, tpr)
        #plt.show()

if __name__ == '__main__':
    print(sys.argv)
    n_extra = int(sys.argv[1])
    lvl = int(sys.argv[2])
    feature_set_dir = sys.argv[3]
    list_of_features = sys.argv[4]
    infile_ped = sys.argv[5]
    stdize = bool(int(sys.argv[6]))
    mtd = sys.argv[7]
    trn_tst_splt = float(sys.argv[8])
    trshold = float(sys.argv[9])
    smote_type = sys.argv[10]
    model_dir = sys.argv[11]
    known_vars = sys.argv[12]
    extra_vars = sys.argv[13]
    myped = ped.Ped(infile_ped)
    myped.addTestFile(field='ind_id', file_pat=os.path.join(feature_set_dir, '%s'))
    myped.ped.dropna(subset=['test'], inplace=True)
    myped.ped.reset_index(inplace=True)
    print('ped shape:')
    print(myped.ped.shape)
    trn = train.TrainTest(known_vars,
                          list_of_features,
                          ['status'],
                          ['descr'])
    trn.stdize = stdize
    trn.threshold = trshold
    trn.readFeatureList()
    trn.readDataSet()
    print('data_set shape is %s' % ' '.join(map(str, trn.data_set.shape)))
    if n_extra > 0:
        trn.readExtraVars(extra_vars, n_extra=n_extra)
    trn.addLabels(level=lvl)
    trn.dropNA('label')
    print('data_set shape is %s' % ' '.join(map(str, trn.data_set.shape)))
    print('label balance is ')
    print(trn.data_set.label.value_counts())
    trn.splitTrainTest(trn_size=trn_tst_splt, over_sample=smote_type)
    print('train_set shape is %s' % ' '.join(map(str, trn.train_set_X.shape)))
    print('test_set shape is %s' % ' '.join(map(str, trn.test_set_X.shape)))
    # set test set equal to train set for final eval
    if trn_tst_splt > .9:
        print('test set is equan to train set')
        trn.test_set_X = trn.train_set_X
        trn.test_set_y = trn.train_set_y
    trn.method = mtd + '_lvl' + str(lvl) +\
                 '_std' + str(trn.stdize) +\
                 '_cut' + str(trn.threshold) +\
                 '_splt' + str(trn_tst_splt) +\
                 '_' + str(n_extra) +\
                 '_' + str(smote_type)
    print('mtd is %s' % mtd)
    if mtd == 'GBM':
 #       n_estiators = 1000
 #       max_depth = 3
 #       learning_rate = 0.01
        n_estimators = 2000
        max_depth = 1
        learning_rate = 0.075
        trn.method += '_' + '_'.join(map(str,
                                   [n_estimators, max_depth, learning_rate]))
        trn.model = GradientBoostingClassifier(n_estimators=n_estimators,
                                               max_depth=max_depth,
                                               learning_rate=learning_rate)
        trn.fitClassifier()
        trn.predictClass(trn.threshold)
    elif mtd == 'LogReg':
        trn.model = LogisticRegression()
        trn.fitClassifier()
        trn.predictClass(trn.threshold)
    elif mtd == 'RF':
        n_estimators = 2000
        max_depth = 1
        trn.method += '_' + '_'.join(map(str, [n_estimators, max_depth]))
        trn.model = RandomForestClassifier(max_depth=max_depth,
                                           n_estimators=n_estimators,
                                           random_state=0)
        trn.fitClassifier()
        trn.predictClass(trn.threshold)
    elif mtd == 'SVM':
        #kernel = 'rbf'
        #C = 1.25
        #class_weight = 'balanced'
        kernel = 'linear'
        C = 0.25
        class_weight = 'balanced'
        trn.method += '_' + '_'.join(map(str,
                                   [kernel, C, class_weight]))
        trn.model = svm.SVC(kernel=kernel, C=C, class_weight=class_weight,
                            probability=True)
#        trn.model = svm.SVC(kernel='linear', C=0.1, class_weight='balanced',
#                            probability=True)
        trn.fitClassifier()
        trn.predictClass(trn.threshold)
    elif mtd == 'OneClassSVM':
        trn.model = svm.OneClassSVM(nu=0.001, kernel='rbf', gamma=0.1,
                                    random_state=0)
        trn.fitClassifier()
        trn.predictClassOneClass()
    elif mtd == 'dl_seq':
        trn.model = Sequential()
        trn.model.add(Dense(256, input_dim=72, init='uniform',
                            activation='tanh'))
#        trn.model.add(Dropout(0.5))
 #       trn.model.add(Dense(256, activation='relu'))
 #       trn.model.add(Dropout(0.5))
        trn.model.add(Dense(1, activation='sigmoid'))
        trn.model.compile(loss='binary_crossentropy', optimizer='sgd')
        trn.model.fit(trn.train_set_X, trn.train_set_y,
                      batch_size=10,
                      nb_epoch=50,
                      shuffle=True,
                      show_accuracy=False)
        x = trn.model.evaluate(trn.test_set_X,
                               trn.test_set_y,
                               batch_size=212,
                               verbose=1,
                               show_accuracy=True)
        trn.pred_y_prob = trn.model.predict_proba(trn.test_set_X)
        trn.pred_y = binarize(trn.pred_y_prob.reshape(1, -1),
                              trn.threshold)[0].astype(int)
        trn.getMetrics()
        trn.saveKerasModel()
        sys.exit(1)
    else:
        sys.exit('Unknown classifier!')
    trn.getMetrics()
    trn.pklModel(model_dir)
