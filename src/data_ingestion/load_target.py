"""
To be written
"""

from pathlib import Path

from utils import read_target


class LoadTarget:
    def __init__(self):
        self.data = None
        self.filePath = None
        self.mode = None

    def load_data(self, path, mode="path"):
        """
        Modes:
        - path: Load AnnData from an H5AD file.
        - object: Receive an existing AnnData object.
        """
        self.mode = mode

        if self.mode == "path":
            self.filePath = Path(path)
            self.data = read_target(self.filePath)

        elif self.mode == "object":
            self.data = path

        else:
            raise ValueError(f"Unsupported loading mode: {self.mode}")

        print(f"Dataset loaded: {self.data.n_obs} cells, {self.data.n_vars} genes")

        return self
