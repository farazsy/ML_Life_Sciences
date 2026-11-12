import pickle
from pathlib import Path

import pandas as pd
import scanpy as sc


def read_text(path):

    path = Path(path)
    extension = path.suffix.lower()

    if extension == ".csv":
        return pd.read_csv(path)

    elif extension in [".xls", ".xlsx"]:
        return pd.read_excel(path)

    else:
        raise ValueError(f"Unsupported file type: {extension}")


def read_embeddings(path):
    """
    Can be cool to investigate and typecast the input I am not sure
    """

    path = Path(path)
    extension = path.suffix.lower()

    if extension not in [".pkl", ".pickle"]:
        raise ValueError(f"Unsupported embedding format: {extension}")

    with open(path, "rb") as file:
        embeddings = pickle.load(file)

    return embeddings


def read_target(path):
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Dataset not found: {path}")

    if path.suffix.lower() != ".h5ad":
        raise ValueError(f"Unsupported target format: {path.suffix}")

    return sc.read_h5ad(path)
