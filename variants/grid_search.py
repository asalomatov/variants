import sys
sys.path.insert(0, '/mnt/xfs1/home/asalomatov/projects/variants/variants')
import ped
import variants
import func
import pandas
import numpy
import os
import features
import train
from multiprocessing import Pool
from sklearn.cross_validation import train_test_split, cross_val_score
from sklearn.ensemble import (RandomForestClassifier,
                              GradientBoostingClassifier,
                              GradientBoostingRegressor)
from sklearn.linear_model import LogisticRegression
from sklearn import svm, metrics
from sklearn.externals import joblib
from sklearn.preprocessing import binarize
from sklearn.preprocessing import scale
from sklearn.decomposition import PCA
from sklearn.grid_search import GridSearchCV, RandomizedSearchCV
from unbalanced_dataset import SMOTE
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
# %matplotlib qt
# sns.set()
# plt.rcParams['font.size'] = 14
# almost_black = '#262626'

# after http://www.codiply.com/blog/hyperparameter-grid-search-across-multiple-models-in-scikit-learn/


class EstimatorSelectionHelper:
    def __init__(self, models, params):
        if not set(models.keys()).issubset(set(params.keys())):
            missing_params = list(set(models.keys()) - set(params.keys()))
            raise ValueError("Some estimators are missing parameters: %s" % missing_params)
        self.models = models
        self.params = params
        self.keys = models.keys()
        self.grid_searches = {}

    def fit(self, X, y, cv=3, n_jobs=1, verbose=1, scoring=None, refit=False):
        for key in self.keys:
            print("Running GridSearchCV for %s." % key)
            model = self.models[key]
            params = self.params[key]
            gs = GridSearchCV(model, params, cv=cv, n_jobs=n_jobs, 
                              verbose=verbose, scoring=scoring, refit=refit)
            gs.fit(X,y)
            self.grid_searches[key] = gs    

    def score_summary(self, sort_by='mean_score'):
        def row(key, scores, params):
            d = {
                 'estimator': key,
                 'min_score': min(scores),
                 'max_score': max(scores),
                 'mean_score': numpy.mean(scores),
                 'std_score': numpy.std(scores),
            }
            return pandas.Series(dict(params.items() + d.items()))

        rows = [row(k, gsc.cv_validation_scores, gsc.parameters) 
                     for k in self.keys
                     for gsc in self.grid_searches[k].grid_scores_]
        df = pandas.concat(rows, axis=1).T.sort([sort_by], ascending=False)

        columns = ['estimator', 'min_score', 'mean_score', 'max_score', 'std_score']
        columns = columns + [c for c in df.columns if c not in columns]

        return df[columns]




if __name__ == '__main__':
    print sys.argv
    n_extra = int(sys.argv[1])
    lvl = int(sys.argv[2])
    feature_set_dir = sys.argv[3]
    list_of_features = sys.argv[4]
    infile_ped = sys.argv[5]
    stdize = bool(int(sys.argv[6]))
    mtd = sys.argv[7]
    trn_tst_splt = float(sys.argv[8])
    trshold = float(sys.argv[9])
    use_class_w = float(sys.argv[10])
    n_jobs = int(sys.argv[11])
    n_folds = int(sys.argv[12])

    known_vars = os.path.join(feature_set_dir, 'fb/known_SNP/fb_known_snp.tsv')
    extra_vars = os.path.join(feature_set_dir, 'fb/all_SNP/fb_all_snp.tsv')
    myped = ped.Ped(infile_ped, ['collection'])
    myped.addTestFile(field='ind_id', file_pat=os.path.join(feature_set_dir, 'fb/all_SNP/%s'))
    myped.ped.dropna(subset=['test'], inplace=True)
    myped.ped.reset_index(inplace=True)
    print myped.ped.shape
    trn = train.TrainTest(known_vars,
                          list_of_features,
                          ['status'],
                          ['descr'])
    trn.stdize = stdize
    trn.threshold = trshold
    trn.readFeatureList()
    trn.readDataSet()
    print 'data_set shape is ', trn.data_set.shape
    if n_extra > 0:
        trn.readExtraVars(extra_vars, n_extra=n_extra)
        print 'data_set shape is ', trn.data_set.shape
    trn.addLabels(level=lvl)
    trn.dropNA('label')
    print 'data_set shape is ', trn.data_set.shape
    trn.splitTrainTest(trn_size=trn_tst_splt)
    print 'train set size ', trn.train_set_X.shape
    print 'testf set size', trn.test_set_X.shape
    trn.method = mtd + '_lvl' + str(lvl) +\
                 '_std' + str(trn.stdize) +\
                 '_cut' + str(trn.threshold) +\
                 '_splt' + str(trn_tst_splt) +\
                 '_' + str(n_extra)
    if mtd == 'GBM':
        param_grid = dict(
            n_estimators=[100, 200, 300, 400, 500, 600, 700, 800, 900, 1000],
            max_depth=[1, 2, 3, 4, 5],
            learning_rate=[0.0001, 0.001, 0.005, 0.01, 0.03, 0.05, 0.075, 0.1, 0.2])
        trn.model = GradientBoostingClassifier(n_estimators=100,
                                               max_depth=1,
                                               learning_rate=0.1)
        grid = GridSearchCV(trn.model, param_grid, cv=n_folds, scoring='roc_auc', n_jobs=n_jobs)
        grid.fit(trn.test_set_X, trn.test_set_y)
        print grid.grid_scores_
        print ' best params , score '
        print grid.best_params_, grid.best_score_
        trn.model = grid.best_estimator_
    elif mtd == 'LogReg':
        trn.model = LogisticRegression()
        trn.fitClassifier()
        trn.predictClass(trn.threshold)
    elif mtd == 'RF':
        trn.model = RandomForestClassifier(max_depth=1,
                                           n_estimators=100,
                                           random_state=0)
        trn.fitClassifier()
        trn.predictClass(trn.threshold)
    elif mtd == 'SVM':
        param_grid = dict(kernel=['linear', 'rbf'],
                          C=[0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5],
                          class_weight=['balanced', None])
        trn.model = svm.SVC(kernel='linear', C=1, probability=True)
        grid = GridSearchCV(trn.model, param_grid, cv=n_folds, scoring='roc_auc', n_jobs=n_jobs)
        grid.fit(trn.test_set_X, trn.test_set_y)
        print grid.grid_scores_
        print ' best params , score '
        print grid.best_params_, grid.best_score_
        trn.model = grid.best_estimator_
    elif mtd == 'OneClassSVM':
        trn.model = svm.OneClassSVM(nu=0.001, kernel='rbf', gamma=0.1,
                                    random_state=0)
        trn.fitClassifier()
        trn.predictClassOneClass()
    else:
        pass
        #sys.exit('Unknown classifier!')
#    trn.getMetrics()


#models1 = { 
#    'RandomForestClassifier': RandomForestClassifier(),
#    'GradientBoostingClassifier': GradientBoostingClassifier(),
#    'SVC': svm.SVC()
#}

models1 = { 
    'GradientBoostingClassifier': GradientBoostingClassifier()
}

params1 = {
    'RandomForestClassifier': dict(n_estimators=[100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]),
    'GradientBoostingClassifier': dict(
            n_estimators=[100, 200, 300, 400, 500, 600, 700, 800, 900, 1000],
            max_depth=[1, 2, 3, 4, 5],
            learning_rate=[0.001, 0.005, 0.01, 0.03, 0.05, 0.075, 0.1, 0.2]),
    'SVC': dict(kernel=['linear', 'rbf'],
                C=[0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5],
                class_weight=['balanced', None])
}

helper1 = EstimatorSelectionHelper(models1, params1)
helper1.fit(trn.train_set_X, trn.train_set_y, cv=n_folds, scoring='roc_auc', n_jobs=n_jobs)

helper1.score_summary(sort_by='mean_score').to_csv('grid_search_GBM.csv',
                                                   index=False)
print helper1.score_summary(sort_by='min_score').head(50)
