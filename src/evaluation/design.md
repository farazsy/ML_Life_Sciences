# Evaluation

## Description

The evaluation module compares predicted gene expression with the
experimental gene expression.

It calculates metrics that measure how well the model predicts the
experimental results.

## Inputs

- y_true: experimental gene expression
- y_pred: predicted gene expression

Both inputs must have the same shape.

## Outputs

The evaluation module returns:

- Mean Squared Error (MSE)
- Root Mean Squared Error (RMSE)
- Pearson correlation

## Metrics

### Mean Squared Error

Measures the average squared difference between predicted and
experimental expression.

Lower values indicate smaller prediction errors.

### Root Mean Squared Error

The square root of MSE.

Lower values indicate smaller prediction errors.

### Pearson Correlation

Measures the linear relationship between predicted and experimental
expression.

Values closer to 1 indicate stronger positive agreement.

## Current Testing

The evaluation module is currently tested using synthetic data.

Once the training and prediction modules are completed, the synthetic
data will be replaced by the actual experimental and predicted
expression data.