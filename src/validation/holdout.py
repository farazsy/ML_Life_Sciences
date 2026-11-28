from src.training.trainer import Trainer
from src.evaluation.evaluate import Evaluate

class HoldoutValidator:
    def __init__(self, trainer: Trainer):
        self.trainer = trainer

    def run(self, X, y, ctrl):
        X_train, X_test, y_train, y_test = self.trainer.split_data(X, y)

        self.trainer.load_model()
        self.trainer.train(X_train, y_train)

        y_pred_train = self.trainer.predict(X_train)
        y_pred_test = self.trainer.predict(X_test)

        evaluator = Evaluate(ctrl)

        return evaluator.report(y_train, y_pred_train, y_test, y_pred_test)