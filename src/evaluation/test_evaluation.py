import numpy as np

from evaluation import Evaluation


# Fake experimental results
y_true = np.array([
    [1, 2, 3],
    [4, 5, 6]
])

# Fake model predictions
y_pred = np.array([
    [1.1, 2.1, 2.9],
    [3.9, 5.2, 5.8]
])


evaluation = Evaluation(y_true, y_pred)

results = evaluation.evaluate()

print("Evaluation results:")
print(results)