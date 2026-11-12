"""
To be written
"""

from pathlib import Path

from utils import read_embeddings, read_text


class LoadTrain:
    def __init__(self):
        self.data = None
        self.filePath = None
        self.mode = None

    def load_data(self, path, mode="embd"):

        self.filePath = Path(path)
        self.mode = mode

        if not self.filePath.exists():
            raise FileNotFoundError(f"Dataset not found: {self.filePath}")

        if self.mode == "text":
            self.data = read_text(self.filePath)

        elif self.mode == "embd":
            self.data = read_embeddings(self.filePath)

        else:
            raise ValueError(f"Invalid read mode: {self.mode}.")

        return self
