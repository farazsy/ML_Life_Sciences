"""
Prepare gene expression targets for the training pipeline.
"""

import scanpy as sc


class BuildTarget:
    def __init__(self, adata):
        self.adata = adata.copy()
        self.y = None

    def filter_cells(self, column="perturbation"):
        """Remove cells without valid perturbation labels."""

        if column not in self.adata.obs.columns:
            raise KeyError(f"Column '{column}' not found in adata.obs")

        mask = self.adata.obs[column].notna() & (self.adata.obs[column] != "*")

        self.adata = self.adata[mask].copy()
        self.adata.obs[column] = self.adata.obs[column].astype(str)

        print(f"Cells after filtering: {self.adata.n_obs}")

        return self

    def format_control(self, column="perturbation"):
        """Standardize control perturbation labels."""

        if column not in self.adata.obs.columns:
            raise KeyError(f"Column '{column}' not found in adata.obs")

        actual_ctrl = ["62(mod)_pBA581", "63(mod)_pBA580", "Gal4-4(mod)_pBA582"]

        self.adata.obs[column] = (
            self.adata.obs[column].replace(actual_ctrl, "ctrl").astype("category")
        )

        return self

    def quality_control(self, min_cells=5, min_counts=200):
        """Filter genes and cells based on QC thresholds."""

        sc.pp.filter_genes(self.adata, min_cells=min_cells)
        sc.pp.filter_cells(self.adata, min_counts=min_counts)

        print(f"After QC: {self.adata.n_obs} cells, {self.adata.n_vars} genes")

        return self

    def normalize(self, mode="target", target_sum=1e6):
        """Normalize total counts per cell."""

        if mode == "target":
            sc.pp.normalize_total(self.adata, target_sum=target_sum)
        elif mode == "default":
            sc.pp.normalize_total(self.adata, target_sum=None)
        else:
            raise ValueError("Mode must be 'target' or 'default'.")

        return self

    def log_transform(self):
        """Apply log1p transformation."""

        sc.pp.log1p(self.adata)

        return self

    def mean_expressions(self, column="perturbation"):
        """Calculate mean gene expression for each perturbation."""

        if column not in self.adata.obs.columns:
            raise KeyError(f"Column '{column}' not found in adata.obs")

        expression_df = self.adata.to_df()

        expression_df[column] = self.adata.obs[column].values

        mean_profiles = expression_df.groupby(column).mean()

        return mean_profiles

    def run_pipeline(
        self,
        min_cells=5,
        min_counts=200,
        normalization="target",
        target_sum=1e6,
        column="perturbation",
    ):
        """Execute the complete target preprocessing pipeline."""

        self.filter_cells(column=column)

        self.format_control(column=column)

        self.quality_control(min_cells=min_cells, min_counts=min_counts)

        self.normalize(mode=normalization, target_sum=target_sum)

        self.log_transform()

        self.y = self.mean_expressions(column=column)

        return self
