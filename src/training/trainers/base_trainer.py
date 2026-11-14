from abc import ABC, abstractmethod


class BaseTrainer(ABC):
    def __init__(self, model_name, model_params=None):
        self.model_name = model_name
        self.model_params = model_params or {}
        self.model = None

    @abstractmethod
    def load_model(self):
        pass

    @abstractmethod
    def train(self, X, y):
        pass

    @abstractmethod
    def predict(self, X):
        pass
