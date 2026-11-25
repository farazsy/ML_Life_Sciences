"""
TBD
"""

from src.data_ingestion.load_target import LoadTarget
from src.data_ingestion.load_train import LoadTrain


class DataIngestion:
    def __init__(self):
        self.train_loader = LoadTrain()
        self.target_loader = LoadTarget()

        self.train_data = None
        self.target_data = None
        self.train_data_type = None

    def load(self, target_path, train_path, train_mode="embd"):

        self.train_data_type = train_mode

        self.train_data = self.train_loader.load_data(train_path, mode=train_mode).data

        self.target_data = self.target_loader.load_data(target_path).data

        return self
