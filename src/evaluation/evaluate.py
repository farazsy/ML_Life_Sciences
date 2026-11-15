import numpy as np
from scipy.stats import pearsonr
from sklearn.metrics import mean_absolute_error, mean_squared_error


class Evaluate:
    def __init__(self, control_mean):
        self.control_mean = np.asarray(control_mean)

    def _evaluate_metrics(self, y_true, y_pred):

        y_true = np.asarray(y_true)
        y_pred = np.asarray(y_pred)

        if y_true.shape != y_pred.shape:
            raise ValueError("Targets and predictions must have matching shapes.")

        true_delta = y_true - self.control_mean
        pred_delta = y_pred - self.control_mean

        rmse_list = []
        mae_list = []
        corr_list = []
        delta_corr_list = []

        for true_row, pred_row, true_d, pred_d in zip(
            y_true, y_pred, true_delta, pred_delta
        ):
            rmse_list.append(np.sqrt(mean_squared_error(true_row, pred_row)))

            mae_list.append(mean_absolute_error(true_row, pred_row))

            if np.std(true_row) > 0 and np.std(pred_row) > 0:
                corr_list.append(pearsonr(true_row, pred_row)[0])

            if np.std(true_d) > 0 and np.std(pred_d) > 0:
                delta_corr_list.append(pearsonr(true_d, pred_d)[0])

        return {
            "rmse": float(np.mean(rmse_list)),
            "mae": float(np.mean(mae_list)),
            "pearson": float(np.mean(corr_list)) if corr_list else np.nan,
            "pearson_delta": float(np.mean(delta_corr_list))
            if delta_corr_list
            else np.nan,
        }

    def report(self, y_train, y_pred_train, y_test, y_pred_test):

        train_metrics = self._evaluate_metrics(y_train, y_pred_train)

        test_metrics = self._evaluate_metrics(y_test, y_pred_test)

        return {"train": train_metrics, "test": test_metrics}
