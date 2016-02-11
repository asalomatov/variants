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
#%matplotlib qt
#sns.set()
#plt.rcParams['font.size'] = 14
#almost_black = '#262626'


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
        self.data_set = self.data_set[~self.data_set.label.isnull()]

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
                                      ['var_id']]

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
        x = x[self.feature_list + [self.y_name]]
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

    def fitClassifier(self):
        # fit estimator
        self.model.fit(self.train_set_X, self.train_set_y)
    #    joblib.dump(self.model, 'estGB_classifier.pkl')
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
        roc_auc = metrics.roc_auc_score(self.test_set_y, self.pred_y_prob)
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
        self.perf_mertics['roc_auc'] = [roc_auc]
        print self.perf_mertics
        #print metrics.classification_report(y_test, y_pred)
        #fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred[:,1]) 
        #plt.plot(fpr, tpr)
        #plt.show()

if __name__ == '__main__':
    import train
#    %matplotlib qt
    sns.set()
    plt.rcParams['font.size'] = 14
    almost_black = '#262626'
    trn = train.Train('/mnt/scratch/asalomatov/data/SSC/wes/feature_sets/fb/known_SNP/fb_known_snp.tsv',
                      '/mnt/xfs1/home/asalomatov/projects/variants/variants/ssc_wes_features.txt',
                      'status')
    trn.readFeatureList()
    trn.readDataSet()
    print 'data_set shape is ', trn.data_set.shape
    trn.readExtraVars('/mnt/scratch/asalomatov/data/SSC/wes/feature_sets/fb/all_SNP/fb_all_snp.tsv', n_extra=10000)
    print 'data_set shape is ', trn.data_set.shape
    trn.addLabels(level=1)
    print 'data_set shape is ', trn.data_set.shape
    trn.data_set.label.isnull().sum()
    print trn.data_set.label.value_counts()
    infile_ped = '/mnt/scratch/asalomatov/data/SSC/SSCped/SSC.ped'
    myped = ped.Ped(infile_ped, ['collection'])
    myped.addTestFile(field='ind_id', file_pat='/mnt/scratch/asalomatov/data/SSC/wes/feature_sets/fb/all_SNP/%s')
    myped.ped.dropna(subset=['test'], inplace=True)
    print myped.ped.shape

    trn.splitTrainTest(trn_size=0.6)

    print '\nGB linear'
    trn.model = GradientBoostingClassifier(n_estimators=100, max_depth=1,
                                           learning_rate=0.1)
    trn.fitClassifier()
    trn.predictClass(threshold=0.5)
    trn.getMetrics()

    print '\nGB depth 3'
    trn.model = GradientBoostingClassifier(n_estimators=100, max_depth=3,
                                           learning_rate=0.1)
    trn.fitClassifier()
    trn.predictClass(threshold=0.5)
    trn.getMetrics()

    print '\nLogisticRegression'
    trn.model = LogisticRegression()
    trn.fitClassifier()
    trn.predictClass(threshold=0.5)
    trn.getMetrics()

    print '\nRandomForest depth 1'
    trn.model = RandomForestClassifier(max_depth=1, n_estimators=100,
                                       random_state=0)
    trn.fitClassifier()
    trn.predictClass(threshold=0.5)
    trn.getMetrics()

    print '\nRandomForest depth 3'
    trn.model = RandomForestClassifier(max_depth=3, n_estimators=100,
                                       random_state=0)
    trn.fitClassifier()
    trn.predictClass(threshold=0.5)
    trn.getMetrics()

    print '\nLinear SVM'
    trn.model = svm.SVC(kernel='linear', C=1, probability=True)
    trn.fitClassifier()
    trn.predictClass(threshold=0.5)
    trn.getMetrics()

    print '\nSVM rbf'
    trn.model = svm.SVC(kernel='rbf', gamma=0.7, C=1, probability=True)
    trn.fitClassifier()
    trn.predictClass(threshold=0.5)
    trn.getMetrics()

    print '\nSVM poly degree 3'
    trn.model = svm.SVC(kernel='poly', degree=3, C=1, probability=True)
    trn.fitClassifier()
    trn.predictClass(threshold=0.5)
    trn.getMetrics()

    sys.exit()


#estGB_regression = GradientBoostingRegressor(n_estimators=2000, max_depth=1, learning_rate=.01)
#len(feature_cols)
#X_train = df_train_set_num[feature_cols].values
#y_train = df_train_set_num[response_col].values
#X_tr, X_te, y_tr, y_te = train_test_split(X_train, y_train)
#X_tr
    sum(estGB_classifier.feature_importances_)




###PCA
    pca = PCA(n_components=3)
    x_vis = pca.fit_transform(Z)
    print pca.explained_variance_ratio_
    x_vis.shape
    x_vis[:5,:]
    x_vis[numpy.array(y) < 0.5, 0].shape
    sum(y < 0.5)
    plot2DimPCA(Z, yZ)
    plot3DimPCA(Z, yZ.astype(int))


    def plot2DimPCA(X, y):
        pca =PCA(n_components = 2)
        x_vis = pca.fit_transform(X)
        print pca.explained_variance_ratio_
        plt.scatter(x_vis[y < 0.5, 0], x_vis[y < 0.5, 1])
        palette = sns.color_palette()
        plt.scatter(x_vis[y < 0.5, 0], x_vis[y < 0.5, 1], label="Class #0", alpha=0.5,
                            edgecolor=almost_black, facecolor=palette[0], linewidth=0.15)
        plt.scatter(x_vis[y == 1, 0], x_vis[y == 1, 1], label="Class #1", alpha=0.5,
                    edgecolor=almost_black, facecolor=palette[2], linewidth=0.15)
        plt.legend()

    def plot3DimPCA(X, y):
        fig = plt.figure(1, figsize=(4, 3))
        plt.clf()
        ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=50)

        plt.cla()
        pca = PCA(n_components=3)
        pca.fit(X)
        X = pca.transform(X)
        print pca.explained_variance_ratio_
        for name, label in [('Negative', 0), ('Positive', 1)]:
            ax.text3D(X[y == label, 0].mean(),                           \
                      X[y == label, 1].mean(),                     \
                      X[y == label, 2].mean(),                     \
                      name,                                              \
                      horizontalalignment='center',                      \
                      bbox=dict(alpha=.5, edgecolor='w', facecolor='w'))
# Reorder the labels to have colors matching the cluster results
        y = numpy.choose(y, [1, 0]).astype(numpy.float)
        ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=y, cmap=plt.cm.spectral)
        x_surf = [X[:, 0].min(), X[:, 0].max(),
                  X[:, 0].min(), X[:, 0].max()]
        y_surf = [X[:, 0].max(), X[:, 0].max(),
                  X[:, 0].min(), X[:, 0].min()]
        x_surf = numpy.array(x_surf)
        y_surf = numpy.array(y_surf)
        v0 = pca.transform(pca.components_[[0]])
        v0 /= v0[-1]
        v1 = pca.transform(pca.components_[[1]])
        v1 /= v1[-1]
        ax.w_xaxis.set_ticklabels([])
        ax.w_yaxis.set_ticklabels([])
        ax.w_zaxis.set_ticklabels([])
        plt.show()

