import numpy as np
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error


class Evaluation:
    """
    Evaluate predicted gene expression against
    experimentally observed gene expression.
    """

    def __init__(self, y_true, y_pred):
        self.y_true = np.asarray(y_true)
        self.y_pred = np.asarray(y_pred)

        if self.y_true.shape != self.y_pred.shape:
            raise ValueError(
                f"Shape mismatch: y_true has shape "
                f"{self.y_true.shape}, but y_pred has shape "
                f"{self.y_pred.shape}"
            )

    def mse(self):
        """Calculate Mean Squared Error."""
        return mean_squared_error(
            self.y_true.flatten(),
            self.y_pred.flatten()
        )

    def rmse(self):
        """Calculate Root Mean Squared Error."""
        return np.sqrt(self.mse())

    def pearson(self):
        """Calculate Pearson correlation."""
        true = self.y_true.flatten()
        pred = self.y_pred.flatten()

        correlation, _ = pearsonr(true, pred)

        return correlation

    def evaluate(self):
        """Calculate all evaluation metrics."""
        return {
            "mse": self.mse(),
            "rmse": self.rmse(),
            "pearson": self.pearson()
        }