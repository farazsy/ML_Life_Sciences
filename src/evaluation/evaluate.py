import numpy as np
from scipy.stats import pearsonr
from sklearn.metrics import mean_absolute_error, mean_squared_error


class Evaluate:
    """
    Evaluate predicted gene expression against experimentally
    observed gene expression, including control-relative metrics.
    """

    def __init__(self, control_mean, metrics=None):
        self.control_mean = np.asarray(control_mean, dtype=float)

        # Default metrics
        self.metrics = metrics or [
            "mse",
            "rmse",
            "mae",
            "pearson",
            "pearson_delta",
        ]

    def _validate(self, y_true, y_pred):

        y_true = np.asarray(y_true, dtype=float)
        y_pred = np.asarray(y_pred, dtype=float)

        if y_true.shape != y_pred.shape:
            raise ValueError(
                f"Shape mismatch: y_true {y_true.shape}, y_pred {y_pred.shape}"
            )

        if y_true.ndim != 2:
            raise ValueError("Targets and predictions must be 2D arrays.")

        # need to check this validation and data shapes
        # if self.control_mean.shape != (y_true.shape[1],):
        # raise ValueError("Control mean must contain one value per target gene.")

        return y_true, y_pred

    def _rowwise_metric(self, y_true, y_pred, metric):

        results = []

        for true_row, pred_row in zip(y_true, y_pred):
            results.append(metric(true_row, pred_row))

        return float(np.mean(results))

    # -------------------------
    # Individual metrics
    # -------------------------

    def mse(self, y_true, y_pred):
        """Mean Squared Error."""
        return self._rowwise_metric(y_true, y_pred, mean_squared_error)

    def rmse(self, y_true, y_pred):
        """Root Mean Squared Error, averaged per perturbation."""
        return self._rowwise_metric(
            y_true, y_pred, lambda true, pred: np.sqrt(mean_squared_error(true, pred))
        )

    def mae(self, y_true, y_pred):
        """Mean Absolute Error."""
        return self._rowwise_metric(y_true, y_pred, mean_absolute_error)

    def pearson(self, y_true, y_pred):
        """Mean Pearson correlation across perturbations."""

        return self._rowwise_metric(y_true, y_pred, self._pearson_row)

    def pearson_delta(self, y_true, y_pred):
        """Pearson correlation of expression changes relative to controls."""

        true_delta = y_true - self.control_mean
        pred_delta = y_pred - self.control_mean

        return self._rowwise_metric(true_delta, pred_delta, self._pearson_row)

    @staticmethod
    def _pearson_row(true, pred):

        if len(true) < 2:
            return np.nan

        if np.std(true) == 0 or np.std(pred) == 0:
            return np.nan

        return float(pearsonr(true, pred)[0])

    # -------------------------
    # Evaluation pipeline
    # -------------------------

    def evaluate(self, y_true, y_pred):

        y_true, y_pred = self._validate(y_true, y_pred)

        results = {}

        for metric_name in self.metrics:
            metric_function = getattr(self, metric_name, None)

            if metric_function is None or not callable(metric_function):
                raise ValueError(f"Unknown metric: {metric_name}")

            results[metric_name] = metric_function(y_true, y_pred)

        return results

    def report(self, y_train, y_pred_train, y_test, y_pred_test):

        train_metrics = self.evaluate(y_train, y_pred_train)
        test_metrics = self.evaluate(y_test, y_pred_test)

        return {"train": train_metrics, "test": test_metrics}
