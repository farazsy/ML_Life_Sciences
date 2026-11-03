"""
Scanpy preprocessing pipeline for perturbation datasets.
"""

from pathlib import Path

import scanpy as sc


class PrePipeline:
    def __init__(self, filePath):
        self.filePath = Path(filePath)
        self.adata = None

    def load_dataset(self):
        """Load AnnData dataset."""

        if not self.filePath.exists():
            raise FileNotFoundError(f"Dataset not found: {self.filePath}")

        self.adata = sc.read_h5ad(self.filePath)

        print(f"Dataset loaded: {self.adata.n_obs} cells, {self.adata.n_vars} genes")

        return self

    def filter_cells(self, column="perturbation"):
        """Remove cells without valid perturbation labels."""

        if column not in self.adata.obs.columns:
            raise KeyError(f"Column '{column}' not found in adata.obs")

        mask = (self.adata.obs[column] != "*") & self.adata.obs[column].notna()

        self.adata = self.adata[mask].copy()

        self.adata.obs["perturbation"] = self.adata.obs["perturbation"].astype(str)

        print(f"Cells after filtering: {self.adata.n_obs}")

        return self

    def format_control(self):
        """Standardize control perturbation labels."""
        actual_ctrl = ["62(mod)_pBA581", "63(mod)_pBA580", "Gal4-4(mod)_pBA582"]

        self.adata.obs["perturbation"] = self.adata.obs["perturbation"].apply(
            lambda x: "ctrl" if x in actual_ctrl else x
        )

        self.adata.obs["perturbation"] = self.adata.obs["perturbation"].astype(
            "category"
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
        else:
            sc.pp.normalize_total(self.adata, target_sum=None)

        return self

    def log_transform(self):
        """Apply log1p transformation."""
        sc.pp.log1p(self.adata)

        return self

    def run_pipeline(
        self, min_cells=5, min_counts=200, normalization="target", target_sum=1e6
    ):
        """Execute the complete preprocessing pipeline."""
        self.load_dataset()

        self.filter_cells()

        self.format_control()

        self.quality_control(min_cells=min_cells, min_counts=min_counts)

        self.normalize(mode=normalization, target_sum=target_sum)

        self.log_transform()

        return self.adata
