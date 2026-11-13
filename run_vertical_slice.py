
import numpy as np
import pandas as pd

from src.data_preparation.align_data import align_data
from src.orchestration.ridge_experiment import run_ridge_experiment
from src.evaluation.evaluate_ridge import evaluate_ridge


def main():
    rng = np.random.default_rng(42)

    #sythetic embeddings
    genes = [f"GENE_{i}" for i in range(20)]
    embeddings = {
        gene: rng.normal(size=8)
        for gene in genes
    }

    #synthetic targets
    expression_genes = [f"EXPR_{i}" for i in range(10)]

    targets = pd.DataFrame(
        rng.normal(size=(20, 10)),
        index=genes,
        columns=expression_genes,
    )

    ##sythetic ctrl
    control_mean = np.zeros(10)

    #align
    X, Y, aligned_genes = align_data(embeddings, targets)

    #train
    result = run_ridge_experiment(X, Y, aligned_genes)

    #eval
    metrics = evaluate_ridge(
        result["actual"],
        result["predictions"],
        control_mean,
    )

    print("\n=== GenePert Vertical Slice ===")
    print(f"Aligned genes: {len(aligned_genes)}")
    print(f"Training genes: {len(result['train_genes'])}")
    print(f"Test genes: {len(result['test_genes'])}")
    print(f"Prediction shape: {result['predictions'].shape}")
    print(f"RMSE: {metrics['rmse']:.4f}")
    print(f"MAE: {metrics['mae']:.4f}")
    print(f"Pearson: {metrics['pers']:.4f}")
    print(f"Pearson delta: {metrics['pers_delta']:.4f}")


if __name__ == "__main__":
    main()
