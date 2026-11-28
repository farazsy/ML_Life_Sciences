"""
TBD
"""

import numpy as np

from src.data_ingestion.data_ingestion import DataIngestion
from src.data_preparation.data_preparation import DataPreparation
from src.validation.validation import Validator
from src.training.trainer import Trainer


class Orch:
    def __init__(self, train_data_path, target_data_path, model="ridge", params=None,
                 validation_method="holdout", n_splits=5, shuffle=True, random_state=42):
        self.train_data_path = train_data_path
        self.target_data_path = target_data_path

        self.X_raw = None
        self.y_raw = None
        self.X_processed = None
        self.y_processed = None

        self.ingestion = DataIngestion()
        self.preparation = DataPreparation()
        self.trainer = Trainer(model, model_params=params, random_state=random_state)
        self.validator = Validator(
            trainer=self.trainer,
            method=validation_method,
            n_splits=n_splits,
            shuffle=shuffle,
            random_state=random_state
        )

        self.report = None

    def run(self, folds=None):
        self.ingestion.load(
            train_path=self.train_data_path, target_path=self.target_data_path
        )

        self.X_raw = self.ingestion.train_data
        self.y_raw = self.ingestion.target_data

        self.preparation.build_dataset(self.X_raw, self.y_raw)

        self.X_processed = np.stack(self.preparation.dataset["embedding"].values)
        self.y_processed = self.preparation.dataset.drop(
            columns=["pert", "embedding"]
        ).to_numpy(dtype=float)

        ctrl = self.preparation.ctrl.drop(columns=["pert"]).to_numpy(dtype=float)
        ids = self.preparation.dataset["pert"].to_numpy()

        self.report = self.validator.run(self.X_processed, self.y_processed, ctrl, ids, folds=folds)

        print(self.report)

        return self.report
