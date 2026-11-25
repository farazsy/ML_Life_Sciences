import numpy as np
import pandas as pd


def align_data(embeddings, targets):
    """
    align pretrained gene embeddings with expression targets.

    embeddings : dict
    targets : pd.DataFrame

    Returns:
    X : np.ndarray
        embeddings
    Y : np.ndarray
        targets
    genes : list
        gene symbols in matching row order
    """
    if not targets.index.is_unique:
        raise ValueError("target gene names must be unique")


    genes = [
        gene for gene in targets.index
        if gene in embeddings
    ]

    if not genes:
        raise ValueError("no matching genes found")

    missing = targets.index.difference(genes).tolist()

    X = np.stack([
        np.asarray(embeddings[gene], dtype=np.float32)
        / np.linalg.norm(embeddings[gene])
        for gene in genes
    ])

    Y = targets.loc[genes].to_numpy(dtype=np.float32)

    print(f"Matched genes: {len(genes)}")
    print(f"Missing embeddings: {len(missing)}")
    print(f"X shape: {X.shape}")
    print(f"Y shape: {Y.shape}")

    return X, Y, genes
