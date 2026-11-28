from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from src.training.trainers.base_model import BaseModel


class PCARegression(BaseModel):
    def __init__(self, n_components=10, alpha=1.0):

        self.n_components = n_components
        self.alpha = alpha
        self.parameters = {"n_components": n_components, "alpha": alpha}

        self.pca = PCA(n_components=n_components)
        self.ridge = Ridge(alpha=alpha)

    def fit(self, X, y):

        # Reduce target dimensionality
        y_reduced = self.pca.fit_transform(y)

        # Train Ridge on the reduced targets
        self.ridge.fit(X, y_reduced)

        return self

    def predict(self, X):

        # Predict reduced targets
        y_reduced = self.ridge.predict(X)

        # Reconstruct original target dimensions
        return self.pca.inverse_transform(y_reduced)
