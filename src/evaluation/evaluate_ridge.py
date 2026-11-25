
import numpy as np
from sklearn.metrics import mean_squared_error, mean_absolute_error
from scipy.stats import pearsonr


def evaluate_ridge(y_true, y_pred, control_mean):
    """
    evaluate predicted expression profiles using RMSE, MAE, pearson correlation and Pearson delta.
    """
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    control_mean = np.asarray(control_mean)

    if y_true.shape != y_pred.shape:
        raise ValueError("predictions and targets must have matching shapes")

    true_delta = y_true - control_mean
    pred_delta = y_pred - control_mean

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
        "pers": float(np.mean(corr_list)) if corr_list else np.nan,
        "pers_delta": float(np.mean(delta_corr_list)) if delta_corr_list else np.nan,
    }
