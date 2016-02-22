import sys
sys.path.insert(0, '/mnt/xfs1/home/asalomatov/projects/variants/variants')
import ped
import variants
import func
import pandas
import numpy
import os
import features
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
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
# %matplotlib qt
# sns.set()
# plt.rcParams['font.size'] = 14
# almost_black = '#262626'


class TrainTest:

    def __init__(self, data_set_file, feature_list_file, y_name):
        self.data_set = None
        self.train_set_X = None
        self.train_set_y = None
        self.train_set_var_id = None
        self.test_set_var_id = None
        self.test_set_X = None
        self.test_set_y = None
        self.pred_y = None
        self.pred_y_prob = None
        self.feature_list = None
        self.feature_list_file = feature_list_file
        self.data_set_file = data_set_file
        self.y_name = y_name
        self.method = None
        self.model = None
        self.perf_mertics = None

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
        if level == 1:
            self.data_set['label'] = None
            self.data_set.ix[self.data_set.status.isin(['Y']), 'label'] = 1
            self.data_set.ix[self.data_set.status.isin(['N']), 'label'] = 0
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
        else:
            sys.exit('Unknown level, exiting...\n')

    def dropNA(self, field):
        self.data_set = self.data_set[~self.data_set[field].isnull()]

    def addVarID(self, df):
        df['var_id'] = df.ind_id +\
                                 '_' +\
                                 df.CHROM.astype(str) +\
                                 '_' +\
                                 df.POS.astype(str)

    def readFeatureList(self):
        """features_file contains names of the columns to be
        used as features, one per line
        """
        with open(self.feature_list_file, 'r') as ff:
            self.feature_list = ff.readlines()
            self.feature_list = [x.strip('\n') for x in self.feature_list]

    def readDataSet(self):
        self.data_set  = pandas.read_csv(self.data_set_file, sep='\t')
        self.addAlleleBalance()
        self.addVarID(self.data_set)
        self.data_set = self.data_set[self.feature_list + [self.y_name] +
                                      ['var_id'] + ['descr']]

    def readTestSet(self, path_to_file):
        #self.data_set_file = path_to_file
        #self.readDataSet()
        #self.test_set_X = self.data_set[self.feature_list].values
        #self.test_set_y = self.data_set.values
        #self.train_set_var_id = list(X_tr.var_id)
        #self.test_set_var_id = list(X_te.var_id)
        pass

    def readExtraVars(self, file_name,  n_extra=50000):
        x  = pandas.read_csv(file_name, sep='\t', nrows=n_extra)
        x = self.addAlleleBalance(x)
        self.addVarID(x)
        x = x[self.feature_list + [self.y_name] + ['var_id'] + ['descr']]
        self.data_set = pandas.concat([self.data_set, x], axis=0)

    def splitTrainTest(self, trn_size=0.5, rnd_state=0):
        X = self.data_set[self.feature_list + ['var_id']]
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

    def data2Test(self):
        c1 = self.data_set.var_id.isin(self.train_set_var_id)
        print sum(c1), 'variants will be removed'
        self.test_set_X = self.data_set[self.feature_list][~c1].values
        self.test_set_y = self.data_set.label[~c1].astype(int).values
        self.test_set_var_id = self.data_set.var_id[~c1].astype(str).values

    def fitClassifier(self):
        # fit estimator
        self.model.fit(self.train_set_X, self.train_set_y)
        model_descr = {'model': self.model,
                       'train_var_id': self.train_set_var_id}
        joblib.dump(model_descr, self.method + '.pkl')
    #    estGB_classifier = joblib.load('estGB_classifier.pkl')
    #    feature_importance = pandas.DataFrame({'contrib': estGB_classifier.feature_importances_ ,'name': feature_columns})
    #    feature_importance.sort(['contrib'], ascending=[False], inplace=True)
    #    feature_importance
    #    #feature_importance.to_excel('feature_contrib.xls', index=False)

    def predictClass(self, threshold=0.5):
        # prediction
        #self.pred_y = self.model.predict(self.test_set_X)
        self.pred_y_prob = self.model.predict_proba(self.test_set_X)[:, 1]
        self.pred_y = binarize(self.pred_y_prob.reshape(1, -1),
                               threshold)[0].astype(int)

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
        if len(numpy.unique(self.test_set_y)) > 1:
            roc_auc = metrics.roc_auc_score(self.test_set_y, self.pred_y_prob)
            self.perf_mertics['roc_auc'] = [roc_auc]
        print self.perf_mertics
        #print metrics.classification_report(y_test, y_pred)
        #fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred[:,1]) 
        #plt.plot(fpr, tpr)
        #plt.show()

if __name__ == '__main__':
    infile_ped = '/mnt/scratch/asalomatov/data/SSC/SSCped/SSC.ped'
    known_vars = '/mnt/scratch/asalomatov/data/SSC/wes/feature_sets/fb/known_SNP/fb_known_snp.tsv'
    list_of_features = '/mnt/xfs1/home/asalomatov/projects/variants/variants/ssc_wes_features.txt'
    myped = ped.Ped(infile_ped, ['collection'])
    myped.addTestFile(field='ind_id', file_pat='/mnt/scratch/asalomatov/data/SSC/wes/feature_sets/fb/all_SNP/%s')
    myped.ped.dropna(subset=['test'], inplace=True)
    myped.ped.reset_index(inplace=True)
    print myped.ped.shape
    n_extra = 10000
    for lvl in [1, 2, 3, 4]:
        trn = TrainTest(known_vars,
                        list_of_features,
                        'status')
        trn.readFeatureList()
        trn.readDataSet()
        print 'data_set shape is ', trn.data_set.shape
        trn.readExtraVars('/mnt/scratch/asalomatov/data/SSC/wes/feature_sets/fb/all_SNP/fb_all_snp.tsv', n_extra=n_extra)
        print 'data_set shape is ', trn.data_set.shape
        trn.addLabels(level=lvl)
        trn.dropNA('label')
        print 'data_set shape is ', trn.data_set.shape
        trn.splitTrainTest(trn_size=0.6)

        print '\nGB additive'
        trn.method = 'GBM_depth1_lvl' + str(lvl) + '_' + str(n_extra)
        trn.model = GradientBoostingClassifier(n_estimators=100, max_depth=1,
                                               learning_rate=0.1)
        trn.fitClassifier()
        trn.predictClass(threshold=0.5)
        trn.getMetrics()

        print '\nLogisticRegression'
        trn.method = 'LogReg_' + str(n_extra)
        trn.model = LogisticRegression()
        trn.fitClassifier()
        trn.predictClass(threshold=0.5)
        trn.getMetrics()

        print '\nRandomForest depth 1'
        trn.method = 'RandForest_depth1_lvl' + str(lvl) + '_' + str(n_extra)
        trn.model = RandomForestClassifier(max_depth=1, n_estimators=100,
                                           random_state=0)
        trn.fitClassifier()
        trn.predictClass(threshold=0.5)
        trn.getMetrics()

        print '\nLinear SVM'
        trn.method = 'LinearSVM_C1_lvl' + str(lvl) + '_' + str(n_extra)
        trn.model = svm.SVC(kernel='linear', C=1, probability=True)
        trn.fitClassifier()
        trn.predictClass(threshold=0.5)
        trn.getMetrics()

        print '\nSVM rbf'
        trn.method = 'rbfSVM_g07_C1_lvl' + str(lvl) + '_' + str(n_extra)
        trn.model = svm.SVC(kernel='rbf', gamma=0.7, C=1, probability=True)
        trn.fitClassifier()
        trn.predictClass(threshold=0.5)
        trn.getMetrics()

        print '\nSVM poly degree 3'
        trn.method = 'polySVM_deg3_C1_lvl' + str(lvl) + '_' + str(n_extra)
        trn.model = svm.SVC(kernel='poly', degree=3, C=1, probability=True)
        trn.fitClassifier()
        trn.predictClass(threshold=0.5)
        trn.getMetrics()
