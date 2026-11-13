import numpy as np
from sklearn.linear_model import Ridge
from sklearn.model_selection import train_test_split


def run_ridge_experiment(
    X,
    Y,
    genes,
    alpha=1.0,
    test_size=0.2,
    random_state=42,
):
    """
    run a minimal GenePert Ridge regression experiment.

    X: embeddings
    Y: targets
    genes: gene names corresponding to rows of X and Y

    Returns: predictions, actual targets, and held-out genes
    """
    X = np.asarray(X)
    Y = np.asarray(Y)

    if len(X) != len(Y) or len(X) != len(genes):
        raise ValueError("X, Y and genes must have matching rows")

    indices = np.arange(len(genes))

    train_idx, test_idx = train_test_split(
        indices,
        test_size=test_size,
        random_state=random_state,
    )

    model = Ridge(alpha=alpha)

    model.fit(X[train_idx], Y[train_idx])
    predictions = model.predict(X[test_idx])

    return {
        "model": model,
        "predictions": predictions,
        "actual": Y[test_idx],
        "train_genes": [genes[i] for i in train_idx],
        "test_genes": [genes[i] for i in test_idx],
    }
