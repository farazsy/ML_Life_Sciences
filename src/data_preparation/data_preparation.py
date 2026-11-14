import pandas as pd

from src.data_preparation.build_target import BuildTarget
from src.data_preparation.build_train import BuildTrain


class DataPreparation:
    def __init__(self, train_mode="load", train_model=None):
        self.train_build = BuildTrain(mode=train_mode, model=train_model)
        self.target_build = BuildTarget()

        self.train_data = None
        self.target_data = None
        self.ctrl = None
        self.dataset = None

    def prepare_train(self, data):
        self.train_build.build(data)
        self.train_data = self.train_build.embeddings

    def prepare_target(self, data):
        self.target_build.run_pipeline(data)
        self.target_data = self.target_build.Y

        df = self.target_data.reset_index()
        self.ctrl = df[df["pert"] == "ctrl"].reset_index(drop=True)

    def build_dataset(self, train_data, target_data):

        # Prepare X and Y
        self.prepare_train(train_data)
        self.prepare_target(target_data)

        # Match perturbation genes with their embeddings
        common = sorted(set(self.train_data.keys()) & set(self.target_data.index))

        # Build X
        X = pd.DataFrame(
            {"pert": common, "embedding": [self.train_data[pert] for pert in common]}
        )

        # Build Y
        Y = self.target_data.loc[common].reset_index(drop=True)

        # Combine X and Y
        self.dataset = pd.concat([X, Y], axis=1)

        return self.dataset
