from sklearn.model_selection import train_test_split

from src.training.models import ModelFactory


class Trainer:
    def __init__(self, model_name, model_params=None, test_size=0.2):

        self.model_name = model_name
        self.model_params = model_params or {}

        self.test_size = test_size

        self.model = None

    def split_data(self, X, y):
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=self.test_size
        )

        return X_train, X_test, y_train, y_test

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
