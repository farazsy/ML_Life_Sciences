import numpy as np
from sklearn.model_selection import KFold

from src.training.trainer import Trainer
from src.evaluation.evaluate import Evaluate

class CrossValidator:
    def __init__(self, trainer: Trainer, n_splits=5, shuffle=True, random_state=42):
        self.trainer = trainer
        self.n_splits = n_splits
        self.shuffle = shuffle
        self.random_state = random_state

    def split_folds(self, X):
        kfold = KFold(n_splits=self.n_splits, shuffle=self.shuffle, 
                      random_state=self.random_state if self.shuffle else None)
        
        return kfold.split(X)
    
    def run_fold(self, X, y, ctrl, train_idx, test_idx):
        X_train = X[train_idx]
        X_test = X[test_idx]

        y_train = y[train_idx]
        y_test = y[test_idx]

        self.trainer.load_model()
        self.trainer.train(X_train, y_train)

        y_pred_train = self.trainer.predict(X_train)
        y_pred_test = self.trainer.predict(X_test)

        evaluator = Evaluate(ctrl)

        return evaluator.report(y_train, y_pred_train, y_test, y_pred_test)
    
    def aggregate_results(self, fold_reports):
        summary = {}

        for split in ["train", "test"]:
            summary[split] = {}

            for metric in fold_reports[0][split]:
                values = np.array([report[split][metric] for report in fold_reports])

                valid_values = values[np.isfinite(values)]

                summary[split][metric] = {
                    "mean": float(np.mean(valid_values)),
                    "std": float(np.std(valid_values, ddof=1)),
                    "n_valid": len(valid_values)
                }

        return summary
    
    def run(self, X, y, ctrl, ids, folds=None):
        fold_reports = []
        fold_details = []
        ids = np.asarray(ids)

        provided_folds = folds is not None

        if folds is None:
            folds = self.split_folds(X)

        for fold, (train_idx, test_idx) in enumerate(folds, start=1):
            report = self.run_fold(X, y, ctrl, train_idx, test_idx)

            fold_reports.append(report)

            fold_details.append({
                "fold": fold,
                "train_ids": ids[train_idx].tolist(),
                "test_ids": ids[test_idx].tolist(),
                "report": report
            })

        return {
            "mode": "cross_validation",
            "n_splits": len(fold_details),
            "fold_source": "provided" if provided_folds else "generated",
            "shuffle": None if provided_folds else self.shuffle,
            "random_state": (None if provided_folds or not self.shuffle
                             else self.random_state),
            "folds": fold_details,
            "summary": self.aggregate_results(fold_reports)
        }