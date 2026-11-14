"""
TBD
"""

import scanpy as sc


class BuildTarget:
    def __init__(self):
        self.data = None
        self.Y = None

    def filter_cells(self, column="perturbation"):
        """Remove cells without valid perturbation labels."""

        if column not in self.data.obs.columns:
            raise KeyError(f"Column '{column}' not found in data.obs")

        mask = (self.data.obs[column] != "*") & self.data.obs[column].notna()

        self.data = self.data[mask].copy()

        self.data.obs[column] = self.data.obs[column].astype(str)

        print(f"Cells after filtering: {self.data.n_obs}")

        return self

    def format_control(self):
        """Standardize control perturbation labels."""
        actual_ctrl = ["62(mod)_pBA581", "63(mod)_pBA580", "Gal4-4(mod)_pBA582"]

        self.data.obs["perturbation"] = self.data.obs["perturbation"].apply(
            lambda x: "ctrl" if x in actual_ctrl else x
        )

        self.data.obs["perturbation"] = self.data.obs["perturbation"].astype("category")

        return self

    def add_pert_column(self):
        """Extract perturbed gene from perturbation label."""

        self.data.obs["pert"] = (
            self.data.obs["perturbation"].astype(str).str.split("_").str[0]
        )

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
        self, data, min_cells=5, min_counts=200, normalization="target", target_sum=1e6
    ):
        """Execute the complete preprocessing pipeline."""

        self.data = data

        self.filter_cells()

        self.format_control()

        self.add_pert_column()

        self.quality_control(min_cells=min_cells, min_counts=min_counts)

        self.normalize(mode=normalization, target_sum=target_sum)

        self.log_transform()

        self.Y = self.mean_expressions()

        return self
