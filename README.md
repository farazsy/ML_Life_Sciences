# GenePert – Machine Learning in Life Sciences

**Team A3**

- Syed Faraz Ahmed
- Vu Ho Tran Anh
- Gwendolyn Elizabeth Howes
- Mohamed Ibrahim Mohamed
- Nilüfer Cagla Yigit

## Project Overview

This project is a modular software implementation and restructuring of **GenePert**, a machine learning approach for predicting gene expression responses to genetic perturbations.

The goal of the project is to organise the original GenePert workflow into a clearer and more reusable software pipeline. The pipeline takes gene information and experimental data as input, prepares the data, trains a prediction model, generates predictions, and evaluates the predictions against experimental results.

The project is being developed as a modular pipeline so that individual components can be developed, tested, and extended independently.

---

## Overall Pipeline

The planned workflow is:

```text
Input Data
    │
    ▼
Data Ingestion
    │
    ▼
Data Preparation
    │
    ├── Load Gene Embeddings
    ├── Build X (Gene Embeddings)
    └── Build Y (Mean Experimental Expression)
    │
    ▼
ML Dataset
    │
    ├── X = Gene Embeddings
    └── Y = Mean Experimental Gene Expression
    │
    ▼
Training
    │
    └── Train GenePert Model
    │
    ▼
Prediction
    │
    └── Generate Predicted Gene Expression
    │
    ▼
Evaluation
    │
    ├── MSE
    ├── RMSE
    └── Pearson Correlation
    │
    ▼
Results
```

An **orchestration/pipeline layer** connects these components and allows the complete workflow to be executed together.

---

## Project Structure

The main implementation is organised under `src/`:

```text
src/
├── data_ingestion/
│   ├── data_ingestion.py
│   ├── load_target.py
│   ├── load_train.py
│   ├── utils.py
│   └── design.md
│
├── data_preparation/
│   ├── build_target.py
│   ├── build_target_1.py
│   ├── build_train.py
│   └── design.md
│
├── training/
│   ├── design.md
│   └── first_draft.py
│
├── evaluation/
│   ├── evaluation.py
│   ├── test_evaluation.py
│   ├── design.md
│   └── first_draft.py
│
└── orchestration/
    ├── design.md
    └── first_draft.py
```

The exact structure may evolve as the individual modules are developed and integrated.

---

## Modules

### 1. Data Ingestion

The data ingestion module is responsible for loading the required input data into the pipeline.

This includes:

- Experimental gene expression data
- Gene information
- Gene embeddings
- Supporting input files

The goal is to keep data loading separate from the later processing and machine learning steps.

---

### 2. Data Preparation

The data preparation stage converts the raw input data into the format required by the machine learning model.

The main components are:

#### Gene Embeddings

Gene embeddings are used to represent genes as numerical feature vectors.

The project is intended to use available published GenePT/GenePert embeddings rather than training a new embedding model as part of this project.

#### Build X

`X` represents the input features used by the machine learning model.

In our pipeline:

```text
X = Gene Embeddings
```

#### Build Y

`Y` represents the experimentally observed gene expression response.

The experimental expression data is grouped by perturbation and mean expression values are calculated to create the target matrix.

Conceptually:

```text
Y = Mean Experimental Gene Expression
```

The resulting machine learning dataset therefore consists of:

```text
X → Gene embeddings
Y → Mean experimental expression
```

---

## 3. Training

The training module is responsible for training the GenePert prediction model using the prepared `X` and `Y` data.

The current approach follows the original GenePert methodology and uses **ridge regression** as the prediction model.

The training stage produces a trained model that can then be used to generate predicted gene expression responses.

---

## 4. Prediction

After training, the model is used to generate predicted gene expression values.

The prediction stage provides the outputs required by the evaluation module.

Conceptually:

```text
Gene embeddings
       ↓
Trained model
       ↓
Predicted gene expression
```

These predictions are then compared against the experimentally observed values.

---

## 5. Evaluation

The evaluation module compares the model predictions with the experimental gene expression values.

The current implementation calculates:

### Mean Squared Error (MSE)

Measures the average squared difference between the predicted and experimental values.

**Lower values indicate smaller prediction errors.**

### Root Mean Squared Error (RMSE)

The square root of MSE.

**Lower values indicate smaller prediction errors.**

### Pearson Correlation

Measures the linear relationship between predicted and experimental values.

**Values closer to 1 indicate stronger positive agreement.**

The evaluation component has been implemented as a reusable `Evaluation` class and currently includes tests using synthetic data.

---

## 6. Orchestration

The orchestration layer connects the individual modules and controls the order in which the pipeline is executed.

The intended flow is:

```text
Data Ingestion
      ↓
Data Preparation
      ↓
Build X / Build Y
      ↓
Training
      ↓
Prediction
      ↓
Evaluation
      ↓
Results
```

This allows the individual components to remain modular while still supporting a complete end-to-end run.

---

## Dataset

The primary dataset used in the project is:

**Adamson et al. (2016)**

Example input:

```text
AdamsonWeissman2016_GSM2406681_10X010.h5ad
```

The data is loaded using the `AnnData`/Scanpy ecosystem.

The preprocessing pipeline includes handling missing or unassigned perturbation labels, quality-control filtering, normalisation, and log transformation before constructing the machine learning targets.

---

## Software Design

The project follows a modular design where each stage has a specific responsibility.

This allows us to:

- Develop components independently
- Test individual components before full integration
- Replace or extend individual stages
- Keep data preparation separate from model training
- Keep evaluation independent from the model implementation
- Run the complete workflow through the orchestration layer

Each major module contains a `design.md` file documenting its planned inputs, outputs, and responsibilities.

---

## Testing

Individual components are tested independently before being connected to the full pipeline.

For example, the evaluation module is currently tested using synthetic experimental and predicted expression data.

This allows us to verify the implementation before the complete training and prediction pipeline is available.

---

## Development Workflow

The team uses Git and GitHub for collaborative development.

Major components can be developed in separate branches before being integrated into the main project.

The project is being developed incrementally:

```text
Design
  ↓
Individual Module
  ↓
Testing
  ↓
Integration
  ↓
End-to-End Pipeline
```

---

## Current Status

The project has progressed from the original GenePert implementation towards a modular pipeline structure.

Current work includes:

- Data ingestion and preparation modules
- Gene embedding handling
- Construction of training and target data
- GenePert model training
- Prediction pipeline
- Evaluation metrics
- Pipeline orchestration
- Independent testing of individual components

The next stage is to continue integrating the individual components and run the complete pipeline using the real experimental dataset.
