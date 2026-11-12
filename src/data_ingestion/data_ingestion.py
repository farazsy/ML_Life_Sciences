"""
TBD
"""

from load_target import LoadTarget
from load_train import LoadTrain


class DataIngestion:
    def __init__(self):
        self.target_loader = LoadTarget()
        self.train_loader = LoadTrain()

        self.adata = None
        self.train_data = None
        self.train_data_type = None

    def load(self, target_path, train_path, train_mode="text"):

        self.adata = self.target_loader.load_data(target_path).adata

        self.train_data = self.train_loader.load_data(train_path, mode=train_mode).data

        self.train_data_type = train_mode

        return self
