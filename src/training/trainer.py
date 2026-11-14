from src.training.models import ModelFactory


class Trainer:
    def __init__(self, model_name, model_params=None):

        self.model_name = model_name
        self.model_params = model_params or {}

        self.model = None

    def load_model(self):

        self.model = ModelFactory.create_model(self.model_name, **self.model_params)

        return self.model

    def train(self, X, y):

        if self.model is None:
            self.load_model()

        self.model.fit(X, y)

        return self

    def predict(self, X):

        if self.model is None:
            raise RuntimeError("Model has not been trained yet.")

        return self.model.predict(X)
