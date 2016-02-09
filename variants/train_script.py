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
from sklearn.ensemble import GradientBoostingClassifier, GradientBoostingRegressor
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from sklearn.preprocessing import binarize
from sklearn.preprocessing import scale
from sklearn.decomposition import PCA
%matplotlib qt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
sns.set()
plt.rcParams['font.size'] = 14
almost_black = '#262626'

def addFeatures(mydf, posND=1):
    mydf['label'] = 0.
    mydf.ix[mydf.status.isin(['Y']), 'label'] = 1.
    if posND == 1: mydf.ix[mydf.status.isin(['ND']), 'label'] = 1.
    mydf.ix[mydf.descr.isin(['both']), 'label'] = 1.
    mydf['allele_balance_offspring'] = mydf['REF_count_offspring']/mydf['DP_offspring']
    mydf['allele_balance_father'] = mydf['REF_count_father']/mydf['DP_father']
    mydf['allele_balance_mother'] = mydf['REF_count_mother']/mydf['DP_mother']
    return mydf

#ssc_snp = pandas.read_csv('/mnt/scratch/asalomatov/data/SSC/wes/feature_sets/fb/all_SNP/fb_all_snp.tsv', \
#       sep='\t')
ssc_snp_all = pandas.read_csv('/mnt/scratch/asalomatov/data/SSC/wes/feature_sets/fb/all_SNP/fb_all_snp.tsv', \
       sep='\t', nrows=50000)
ssc_snp = pandas.read_csv('/mnt/scratch/asalomatov/data/SSC/wes/feature_sets/fb/known_SNP/fb_known_snp.tsv', \
       sep='\t')
ssc_snp.shape
ssc_snp.status.value_counts()
ssc_snp.descr.value_counts()
ssc_snp.dtypes
ssc_snp = addFeatures(ssc_snp)
ssc_snp_all = addFeatures(ssc_snp_all)
ssc_snp.label.value_counts()

feature_columns = [i for i in ssc_snp.columns if i not in ['CHROM', 'POS', 'REF_base_offspring', 'ALT_base_offspring',\
        'REF_base_father', 'ALT_base_father', 
        'REF_base_mother', 'ALT_base_mother', 
        'status', 'label', 'ind_id', 'descr']]
len(feature_columns)
X = ssc_snp[feature_columns].values
y = ssc_snp.label.values.astype(int)
XX = ssc_snp_all[feature_columns].values
yy = ssc_snp_all.label.values
Z = pandas.concat([ssc_snp[feature_columns], ssc_snp_all[feature_columns]], axis=0).values
yZ = numpy.concatenate([y, yy])
yZ - yZ.astype(int)
Z.shape
yZ.shape
y.sum()
(y == 0).sum()
yy.sum()
type(y)
X[ y < .5]

X_train, X_test, y_train, y_test = train_test_split(Z, yZ, train_size=.8, random_state=0)
X_train.shape
X_test.shape
X_test.columns
sum(y_train)
sum(y_test)
y_test.mean()

# fit estimator
estGB_classifier = GradientBoostingClassifier(n_estimators=100, max_depth=1, learning_rate=.1)
#estGB_regression = GradientBoostingRegressor(n_estimators=2000, max_depth=1, learning_rate=.01)
#len(feature_cols)
#X_train = df_train_set_num[feature_cols].values
#y_train = df_train_set_num[response_col].values
#X_tr, X_te, y_tr, y_te = train_test_split(X_train, y_train)
#X_tr
estGB_classifier.fit(X_train, y_train)
feature_importance = pandas.DataFrame({'contrib': estGB_classifier.feature_importances_ ,'name': feature_columns})
feature_importance.sort(['contrib'], ascending=[False], inplace=True)
feature_importance
#feature_importance.to_excel('feature_contrib.xls', index=False)
sum(estGB_classifier.feature_importances_)

# prediction
y_pred = estGB_classifier.predict(X_test)
y_pred_prob = estGB_classifier.predict_proba(X_test)
y_pred.shape
y_pred[:5,:]
metrics.accuracy_score(y_test, y_pred)  # accuracy
max(y_test.mean(), 1 - y_test.mean())  # null accuracy
metrics.confusion_matrix(y_test, y_pred)
print metrics.classification_report(y_test, y_pred)

#fpr, tpr, thresholds = metrics.roc_curve(y_test, binarize(y_pred, .7)) 
fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred[:,1]) 
plt.plot(fpr, tpr)



### PCA
pca =PCA(n_components = 3)
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
    plt.scatter(x_vis[y==1, 0], x_vis[y==1, 1], label="Class #1", alpha=0.5, 
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




