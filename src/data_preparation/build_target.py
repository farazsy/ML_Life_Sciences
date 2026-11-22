"""
TBD
"""

import numpy as np
import scanpy as sc

from src.data_preparation.utils import calculate_mad


class BuildTarget:
    def __init__(self):
        self.data = None
        self.Y = None

    def set_data(self, data, copy=False):
        """Set the AnnData object used by the pipeline."""

        if data is None:
            raise ValueError("Input data cannot be None")

        if not hasattr(data, "obs") or not hasattr(data, "X"):
            raise TypeError("Input data must be an AnnData object")

        self.data = data.copy() if copy else data

        self.data.var_names_make_unique()

        self.Y = None

        return self

    def filter_cells(self, column="perturbation", drop_labels=("*",)):
        """Remove cells without valid perturbation labels."""

        if column not in self.data.obs.columns:
            raise KeyError(f"Column '{column}' not found in data.obs")

        labels = self.data.obs[column].astype("string").str.strip()

        mask = labels.notna() & labels.ne("") & ~labels.isin(drop_labels)

        self.data = self.data[mask].copy()
        self.data.obs[column] = labels.loc[mask].astype(str)

        print(f"Cells after filtering: {self.data.n_obs}")

        return self

    def format_control(
        self, column="perturbation", control_labels=None, ctrl_label="ctrl"
    ):
        """Standardize control perturbation labels."""
        if column not in self.data.obs.columns:
            raise KeyError(f"Column '{column}' not found in data.obs")

        if control_labels is None:
            control_labels = ["62(mod)_pBA581", "63(mod)_pBA580", "Gal4-4(mod)_pBA582"]

        self.data.obs[column] = (
            self.data.obs[column]
            .replace(dict.fromkeys(control_labels, ctrl_label))
            .astype("category")
        )

        return self

    def add_pert_column(
        self, column="perturbation", output_column="pert", separator="_"
    ):
        """Extract perturbed gene from perturbation label."""

        if column not in self.data.obs.columns:
            raise KeyError(f"Column '{column}' not found in data.obs")

        labels = self.data.obs[column].astype(str)

        self.data.obs[output_column] = labels.str.split(separator, n=1).str[0]

        return self

    def mad_filter(self, max_mad=5):
        """Remove cells with unusually high total counts using MAD."""

        X = self.data.X

        n_counts = np.asarray(X.sum(axis=1)).ravel()

        median, mad = calculate_mad(n_counts)

        upper_limit = median + max_mad * mad

        mask = n_counts <= upper_limit

        self.data = self.data[mask].copy()

        print(f"Cells removed by MAD: {(~mask).sum()}")
        print(f"Cells remaining: {self.data.n_obs}")

        return self

    def quality_control(self, min_cells=5, min_counts=200):
        """Filter genes and cells based on QC thresholds."""
        sc.pp.filter_genes(self.data, min_cells=min_cells)

        sc.pp.filter_cells(self.data, min_counts=min_counts)

        print(f"After QC: {self.data.n_obs} cells, {self.data.n_vars} genes")

        return self

    def normalize(self, mode="target", target_sum=1e6):
        """Normalize total counts per cell."""
        if mode == "target":
            sc.pp.normalize_total(self.data, target_sum=target_sum)
        else:
            sc.pp.normalize_total(self.data, target_sum=None)

        return self

    def log_transform(self):
        """Apply log1p transformation."""
        sc.pp.log1p(self.data)

        return self

    def mean_expressions(self, column="pert"):
        if column not in self.data.obs.columns:
            raise KeyError(f"Column '{column}' not found in data.obs")

        # Reconsider the numpy approach
        expression_df = self.data.to_df()

        expression_df[column] = self.data.obs[column].values

        mean_profiles = expression_df.groupby(column).mean()

        return mean_profiles

    def run_pipeline(
        self,
        data,
        control_labels=None,
        ctrl_label="ctrl",
        perturbation_column="perturbation",
        gene_column="pert",
        separator="_",
        mad_max=5,
        min_cells=5,
        min_counts=200,
        normalization="target",
        target_sum=1e6,
        copy_data=True,
    ):
        """Execute the complete preprocessing pipeline."""

        self.set_data(data, copy=copy_data)

        # Label processing
        self.filter_cells(column=perturbation_column)

        self.format_control(
            column=perturbation_column,
            control_labels=control_labels,
            ctrl_label=ctrl_label,
        )

        self.add_pert_column(
            column=perturbation_column, output_column=gene_column, separator=separator
        )

        # Quality control
        self.mad_filter(max_mad=mad_max)
        self.quality_control(min_cells=min_cells, min_counts=min_counts)

        # Expression preprocessing
        self.normalize(mode=normalization, target_sum=target_sum)
        self.log_transform()

        # Generate target profiles
        self.Y = self.mean_expressions()

        return self
