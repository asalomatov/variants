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

