# Machine Learning in Life Sciences | Group A3

**FU Berlin | MSc Data Science**

## Team Members
- Syed Faraz Ahmed
- Mohamed ibrahim mohamed Abdelhamid

## Project Overview

This project focuses on applying machine learning methods to perturbation data in life sciences.

Our work includes:
- Preprocessing the Adamson et al. (2016) dataset.
- Implementing and adapting the GenePert approach.
- Using gene embeddings to predict gene expression changes after perturbations.
- Evaluating model performance using appropriate metrics and visualisations.

## Repository Structure

```text
ML_Life_Sciences/
│
├── genepert/
│   ├── GenePertExperiment.py
│   └── utils.py
│
├── preprocessing/
│   ├── __init__.py
│   └── preprocessing.py
│
├── genepert-k562-demo.ipynb
├── test.py
└── .gitignore
```

## Requirements

The project uses Python and the following main libraries:

- Scanpy
- AnnData
- PyTorch
- Scikit-learn
- Pandas
- NumPy
- Matplotlib
- Seaborn
- UMAP

Install the dependencies using:

```bash
pip install scanpy anndata torch scikit-learn pandas numpy matplotlib seaborn umap-learn
```

## Dataset

The project uses the Adamson et al. (2016) perturbation dataset.

The original and processed `.h5ad` files are not included in this repository because of their large file sizes.

Dataset download:

[AdamsonWeissman2016 dataset](https://zenodo.org/record/7041849/files/AdamsonWeissman2016_GSM2406681_10X010.h5ad)

## Running the Preprocessing

After placing the dataset in the project root directory, activate your Python environment and run:

```bash
python test.py
```

This runs the preprocessing pipeline and generates `processed_dataset.h5ad`.

## GenePert

The `genepert/` directory contains the GenePert experiment implementation and utility functions. The demo notebook is also included as a reference for the original approach.

Further adaptations and experiments will be documented as the project progresses.
