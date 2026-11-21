from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge

from src.training.custom_models.base_model import BaseModel


class PCARegression(BaseModel):


    def __init__(self, n_components=10, alpha=1.0):
        self.n_components = n_components
        self.alpha = alpha

    def fit(self, X, y):
        self.pca = PCA(n_components=self.n_components)
        y_reduced = self.pca.fit_transform(y)

        self.ridge = Ridge(alpha=self.alpha)
        self.ridge.fit(X, y_reduced)

        return self

    def predict(self, X):
        y_reduced = self.ridge.predict(X)

        return self.pca.inverse_transform(y_reduced)