from sklearn.linear_model import Ridge
from sklearn.neighbors import KNeighborsRegressor
from src.training.custom_models.pca_regression import PCARegression


class ModelFactory:
    @staticmethod
    def create_model(model_name, **kwargs):

        models = {"ridge": Ridge,
                  "knn": KNeighborsRegressor,
                  "pca_ridge": PCARegression}

        if model_name not in models:
            raise ValueError(f"Unknown model: {model_name}")

        return models[model_name](**kwargs)
