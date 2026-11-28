"""
TBD
"""

import json
import numpy as np

from src.data_ingestion.data_ingestion import DataIngestion
from src.data_preparation.data_preparation import DataPreparation
from src.evaluation.evaluate import Evaluate
from src.training.trainer import Trainer


class Orch:
    def __init__(self, config):
        self.config = config

        self.X_raw = None
        self.y_raw = None
        self.X_processed = None
        self.y_processed = None
        self.X_train = None
        self.X_test = None
        self.y_train = None
        self.y_test = None
        self.y_pred_test = None
        self.y_pred_train = None

        self.ingestion = DataIngestion()
        self.preparation = DataPreparation(**config["DataPreparation"])
        self.trainer = Trainer(**config["Trainer"])
        self.evaluate = None

        self.report = None

    @classmethod
    def from_json(cls, path):
        with open(path) as f:
            config = json.load(f)

        return cls(config)

    def run(self):
        self.ingestion.load(**self.config["DataIngestion"])

        self.X_raw = self.ingestion.train_data
        self.y_raw = self.ingestion.target_data

        self.preparation.build_dataset(self.X_raw, self.y_raw)

        self.X_processed = np.stack(self.preparation.dataset["embedding"].values)
        self.y_processed = self.preparation.dataset.drop(
            columns=["pert", "embedding"]
        ).to_numpy(dtype=float)

        self.X_train, self.X_test, self.y_train, self.y_test = self.trainer.split_data(
            self.X_processed, self.y_processed
        )

        ctrl = self.preparation.ctrl.drop(columns=["pert"]).to_numpy(dtype=float)

        self.trainer.train(self.X_train, self.y_train)

        self.y_pred_train = self.trainer.predict(self.X_train)
        self.y_pred_test = self.trainer.predict(self.X_test)

        self.evaluate = Evaluate(ctrl, **self.config["Evaluate"])

        self.report = self.evaluate.report(
            self.y_train, self.y_pred_train, self.y_test, self.y_pred_test
        )

        print(self.report)
