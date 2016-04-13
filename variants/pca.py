import sys
#sys.path.insert(0, '/mnt/xfs1/home/asalomatov/projects/variants/variants')
import numpy
import train
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
sns.set()
plt.rcParams['font.size'] = 14
almost_black = '#262626'


def plot2DimPCA(X, y, file_name=None):
    pca = PCA(n_components=2)
    x_vis = pca.fit_transform(X)
    print 'pca explained variance ratio :', pca.explained_variance_ratio_
    plt.scatter(x_vis[y < 0.5, 0], x_vis[y < 0.5, 1])
    palette = sns.color_palette()
    plt.scatter(x_vis[y < 0.5, 0], x_vis[y < 0.5, 1], label="Class #0", alpha=0.5,
                        edgecolor=almost_black, facecolor=palette[0], linewidth=0.15)
    plt.scatter(x_vis[y == 1, 0], x_vis[y == 1, 1], label="Class #1", alpha=0.5,
                edgecolor=almost_black, facecolor=palette[2], linewidth=0.15)
    plt.legend()
    if file_name is not None:
        plt.figure().savefig(file_name, bbox_inches='tight')


def plot3DimPCA(X, y, file_name=None):
    fig = plt.figure(1, figsize=(4, 3))
    plt.clf()
    ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=50)

    plt.cla()
    pca = PCA(n_components=3)
    pca.fit(X)
    X = pca.transform(X)
    print pca.explained_variance_ratio_
    for name, label in [('Negative', 0), ('Positive', 1)]:
        ax.text3D(X[y == label, 0].mean(),
                  X[y == label, 1].mean(),
                  X[y == label, 2].mean(),
                  name,
                  horizontalalignment='center',
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
    if file_name is not None:
        plt.figure().savefig(file_name, bbox_inches='tight')


if __name__ == '__main__':
    # run PCA
#    infile_ped = '/mnt/scratch/asalomatov/data/SSC/SSCped/SSC.ped'
    known_vars = '/mnt/scratch/asalomatov/data/SSC/wes/feature_sets_01/fb/known_SNP/fb_known_snp.tsv'
    extra_vars = '/mnt/scratch/asalomatov/data/SSC/wes/feature_sets_01/fb/all_SNP/fb_all_snp.tsv'
    list_of_features = '/mnt/xfs1/home/asalomatov/projects/variants/variants/ssc_wes_features_noINDEL_noDP.txt'
#    myped = ped.Ped(infile_ped, ['collection'])
#    myped.addTestFile(field='ind_id', file_pat='/mnt/scratch/asalomatov/data/SSC/wes/feature_sets_01/fb/all_SNP/%s')
#    myped.ped.dropna(subset=['test'], inplace=True)
#    myped.ped.reset_index(inplace=True)
#    print myped.ped.shape

    n_extra = int(sys.argv[1])
    lvl = int(sys.argv[2])
    dim = int(sys.argv[3])
    trn = train.TrainTest(known_vars,
                          list_of_features,
                          ['status'], ['descr'])
    trn.readFeatureList()
    trn.readDataSet()
    print 'data_set shape is ', trn.data_set.shape
    trn.readExtraVars(extra_vars, n_extra=n_extra)
    print 'data_set shape is ', trn.data_set.shape
    trn.addLabels(level=lvl)
    trn.dropNA('label')
#    sys.exit(1)
    trn.data2Test()
    print 'test_set_X shape is', trn.test_set_X.shape
    if dim ==2:
        plot2DimPCA(trn.test_set_X, trn.test_set_y)
    elif dim ==3:
        plot3DimPCA(trn.test_set_X, trn.test_set_y)
    else:
        sys.exit('\nUnknown dim!')
    plt.show()

