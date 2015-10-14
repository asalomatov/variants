from sklearn import datasets
iris = datasets.load_iris()
diab = datasets.load_diabetes()
print iris.data
print doab.data
type(iris.data)



from sklearn.datasets import make_hastie_10_2
from sklearn.cross_validation import train_test_split
from sklearn.ensemble import GradientBoostingClassifier, GradientBoostingRegressor



# generate synthetic data from ESLII - Example 10.2
X, y = make_hastie_10_2(n_samples=5000)
X_train, X_test, y_train, y_test = train_test_split(X, y)

# fit estimator
est = GradientBoostingClassifier(n_estimators=200, max_depth=3)
est.fit(X_train, y_train)

# predict class labels
pred = est.predict(X_test)

# score on test data (accuracy)
acc = est.score(X_test, y_test)
print('ACC: %.4f' % acc)

# predict class probabilities
est.predict_proba(X_test)[0]


