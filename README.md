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
├── scripts/
│   └── preprocess_adamson.py      # CLI: raw Adamson h5ad -> GenePert-ready h5ad
│
├── src/
│   ├── data_ingestion/            # loaders (h5ad, embeddings, text)
│   └── data_preparation/
│       ├── build_target.py        # generic QC / HVG / normalization / pseudo-bulk
│       └── adamson.py             # Adamson guide parsing, guide QC, symbol fixes
│
├── genepert-k562-demo.ipynb
├── requirements.txt
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
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## Dataset

The project uses the Adamson et al. (2016) perturbation dataset.

The original and processed `.h5ad` files are not included in this repository because of their large file sizes.

Dataset download:

[AdamsonWeissman2016 dataset](https://zenodo.org/record/7041849/files/AdamsonWeissman2016_GSM2406681_10X010.h5ad)

## Running the Preprocessing

After placing the dataset in the project root directory, activate your Python environment and run:

```bash
python scripts/preprocess_adamson.py          # or: python test.py
python scripts/preprocess_adamson.py --help   # all thresholds are configurable
```

Outputs in `data/processed/` (git-ignored):

- `adamson_genepert.h5ad`: log1p(CP10k) expression, 5,000 control-derived HVGs + perturbed genes,
  raw counts in `layers['counts']`, `obs['condition']` as `GENE+ctrl` / `ctrl` (GenePert/GEARS format)
- `guide_qc.csv`: per-guide cell count, knockdown ratio, sibling-guide agreement, kept/flagged
- `adamson_genepert_pseudobulk.csv`: mean expression per condition

Pipeline (rationale in `GenePT_GenePert_study_notes.md`, §5–7):

1. Drop cells with no guide (`NaN`) or ambiguous guide (`*`); pool `63(mod)` + `Gal4-4(mod)` as `ctrl`, drop `62(mod)`.
2. Cell QC: mito > 15 %, < 500 genes, total counts > median + 5·MAD (doublets); genes detected in < 10 cells (targets kept).
3. Guide QC: drop guides whose target is not knocked down (mean target / ctrl > 0.7), or whose
   knockdown can't be measured and whose profile disagrees with sibling guides (r < 0.5).
4. Merge guides into genes; keep genes with ≥ 50 cells.
5. Rename old hg19 symbols to current HGNC (`TARS` → `TARS1`, `C7orf26` → `INTS15`, …; use `--no-harmonize` to keep them) and fix `IER3IP1` by Ensembl ID.
6. HVGs (`seurat_v3`) on control cells only, so they don't leak test perturbations; add back perturbed genes.
7. Normalize to 10,000 counts over all genes, log1p, subset genes. No batch correction or scaling.

Default run: 60,797 cells (7,206 ctrl), 82 perturbations, 5,062 genes.

## GenePert

The `genepert/` directory contains the GenePert experiment implementation and utility functions. The demo notebook is also included as a reference for the original approach.

Further adaptations and experiments will be documented as the project progresses.
