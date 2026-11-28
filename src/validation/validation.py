from src.training.trainer import Trainer
from src.validation.holdout import HoldoutValidator
from src.validation.cross_validation import CrossValidator

class Validator:
    def __init__(self, trainer:Trainer, method="holdout", n_splits=5, shuffle=True, random_state=42):
        self.trainer = trainer
        self.method = method
        self.n_splits = n_splits
        self.shuffle = shuffle
        self.random_state = random_state

    def select_strategy(self):
        if self.method == "holdout":
            return HoldoutValidator(self.trainer)
        
        elif self.method == "cross_validation":
            return CrossValidator(trainer=self.trainer, n_splits=self.n_splits, 
                                  shuffle=self.shuffle, random_state=self.random_state)
        
        else:
            raise ValueError(f"Unknown validation method: {self.method}")
        
    def run(self, X, y, ctrl, ids=None, folds=None):
        strategy = self.select_strategy()

        if self.method == "holdout":
            return strategy.run(X, y, ctrl)
        
        elif self.method == "cross_validation":
            if ids is None:
                raise ValueError("Cross-validation requires perturbation IDs.")
            
            return strategy.run(X, y, ctrl, ids, folds=folds)