from sklearn.model_selection import train_test_split

from src.training.model_factory import ModelFactory


class Trainer:
    def __init__(self, model_name, model_params=None, test_size=0.2, random_state=42):

        self.model_name = model_name
        self.model_params = model_params or {}

        self.test_size = test_size
        self.random_state = random_state

        self.model = None

    def split_data(self, X, y):

        return train_test_split(
            X, y, test_size=self.test_size, random_state=self.random_state
        )

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
