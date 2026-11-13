"""
Prepare gene expression targets (pseudo-bulk profiles per perturbation).

Generic, dataset-agnostic steps. Dataset-specific label handling (e.g. Adamson
guide -> gene mapping) lives in its own module (see ``adamson.py``).

Typical order:
    filter_cells -> format_control -> qc_cells -> quality_control
    -> select_hvgs (on raw counts) -> normalize -> log_transform
    -> subset_genes -> set_condition -> mean_expressions
"""

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


class BuildTarget:
    def __init__(self, adata, copy=False):
        # Copying the raw Adamson matrix doubles memory (~2 GB), so default to in-place.
        self.adata = adata.copy() if copy else adata
        self.Y = None
        self.keep_genes = None
        self._counts = None

    # ------------------------------------------------------------------ labels
    def filter_cells(self, column="perturbation", drop_labels=("*",)):
        """Remove cells without a valid perturbation label (NaN or ambiguous)."""
        self._require(column)

        labels = self.adata.obs[column]
        mask = labels.notna() & ~labels.isin(drop_labels)

        self.adata = self.adata[mask.values].copy()
        self.adata.obs[column] = self.adata.obs[column].astype(str)

        print(f"Cells after label filtering: {self.adata.n_obs}")
        return self

    def format_control(self, control_labels, column="perturbation", ctrl_label="ctrl"):
        """Map all control guide labels to a single control label."""
        self._require(column)

        self.adata.obs[column] = (
            self.adata.obs[column]
            .astype(str)
            .replace(list(control_labels), ctrl_label)
            .astype("category")
        )
        return self

    # ---------------------------------------------------------------------- QC
    def qc_cells(self, max_pct_mito=15.0, min_genes=500, max_counts_mad=5.0):
        """
        Cell-level QC on raw counts.

        - drop cells with percent_mito > max_pct_mito
        - drop cells with fewer than min_genes detected genes
        - drop likely doublets: total counts > median + max_counts_mad * MAD
        """
        X = self.adata.X
        n_counts = np.asarray(X.sum(axis=1)).ravel()
        n_genes = np.asarray((X > 0).sum(axis=1)).ravel()

        if "percent_mito" in self.adata.obs:
            pct_mito = self.adata.obs["percent_mito"].to_numpy()
        else:
            mito = self.adata.var_names.str.upper().str.startswith("MT-")
            pct_mito = 100 * np.asarray(X[:, mito].sum(axis=1)).ravel() / n_counts

        median = np.median(n_counts)
        mad = np.median(np.abs(n_counts - median))
        max_counts = median + max_counts_mad * mad

        mask = (pct_mito <= max_pct_mito) & (n_genes >= min_genes) & (n_counts <= max_counts)

        print(
            f"Cell QC: mito>{max_pct_mito}%: {(pct_mito > max_pct_mito).sum()}, "
            f"genes<{min_genes}: {(n_genes < min_genes).sum()}, "
            f"counts>{max_counts:.0f}: {(n_counts > max_counts).sum()} "
            f"-> removing {(~mask).sum()} cells"
        )

        self.adata = self.adata[mask].copy()
        return self

    def quality_control(self, min_cells=10, min_counts=None, protect_genes=()):
        """
        Drop genes detected in fewer than min_cells cells (and optionally low-count cells).
        Genes in protect_genes (e.g. perturbation targets) are never dropped.
        """
        n_cells = np.asarray((self.adata.X > 0).sum(axis=0)).ravel()
        gene_mask = (n_cells >= min_cells) | self.adata.var_names.isin(list(protect_genes))
        self.adata = self.adata[:, gene_mask].copy()
        if min_counts is not None:
            sc.pp.filter_cells(self.adata, min_counts=min_counts)

        print(f"After QC: {self.adata.n_obs} cells, {self.adata.n_vars} genes")
        return self

    # ----------------------------------------------------------- HVG selection
    def select_hvgs(
        self,
        n_top_genes=5000,
        column="condition",
        ctrl_label="ctrl",
        control_only=True,
        add_genes=(),
        flavor="seurat_v3",
    ):
        """
        Choose HVGs on RAW counts and remember them (plus ``add_genes``) for ``subset_genes``.

        With control_only=True HVGs are computed on control cells only, so the gene set
        does not depend on which perturbations end up in train or test folds.
        Must run before ``normalize`` (seurat_v3 expects counts). The actual subset is
        applied later so that normalization uses total counts over all genes.
        """
        self._require(column)

        if control_only:
            ref = self.adata[self.adata.obs[column] == ctrl_label].copy()
        else:
            ref = self.adata.copy()

        sc.pp.highly_variable_genes(ref, n_top_genes=n_top_genes, flavor=flavor)
        hvgs = ref.var_names[ref.var["highly_variable"]]

        extra = [g for g in add_genes if g in self.adata.var_names and g not in set(hvgs)]
        missing = [g for g in add_genes if g not in self.adata.var_names]
        if missing:
            print(f"Warning: {len(missing)} add_genes not in var: {missing}")

        keep = set(hvgs) | set(extra)
        self.keep_genes = [g for g in self.adata.var_names if g in keep]  # keep var order
        self.adata.var["highly_variable"] = self.adata.var_names.isin(hvgs)

        # Stash raw counts of the kept genes for layers['counts'].
        self._counts = self.adata[:, self.keep_genes].X.copy()

        print(f"HVGs: {len(hvgs)} + {len(extra)} added perturbed genes = {len(self.keep_genes)}")
        return self

    def subset_genes(self, keep_counts_layer=True):
        """Apply the gene subset chosen by select_hvgs."""
        if self.keep_genes is None:
            raise RuntimeError("Call select_hvgs before subset_genes.")

        self.adata = self.adata[:, self.keep_genes].copy()
        if keep_counts_layer and self._counts is not None:
            self.adata.layers["counts"] = self._counts
        self._counts = None
        return self

    # ------------------------------------------------------- normalization
    def normalize(self, target_sum=1e4):
        """Normalize total counts per cell (target_sum=None -> median library size)."""
        sc.pp.normalize_total(self.adata, target_sum=target_sum)
        return self

    def log_transform(self):
        """Apply log1p transformation."""
        sc.pp.log1p(self.adata)
        return self

    # ---------------------------------------------------------- conditions
    def set_condition(self, gene_column="gene", ctrl_label="ctrl", out_column="condition"):
        """Build GEARS/GenePert style labels: 'ctrl' and 'GENE+ctrl'."""
        self._require(gene_column)

        genes = self.adata.obs[gene_column].astype(str)
        self.adata.obs[out_column] = pd.Categorical(
            np.where(genes == ctrl_label, ctrl_label, genes + "+" + ctrl_label)
        )
        return self

    def mean_expressions(self, column="condition"):
        """Pseudo-bulk: mean expression per condition (conditions x genes DataFrame)."""
        self._require(column)
        return group_means(self.adata.X, self.adata.obs[column], self.adata.var_names)

    # ---------------------------------------------------------- pipeline
    def run_pipeline(
        self,
        control_labels=(),
        column="perturbation",
        min_cells=10,
        target_sum=1e4,
        n_top_genes=5000,
    ):
        """Generic pipeline (no dataset-specific guide QC). Uses `column` as condition."""
        self.filter_cells(column=column)
        self.format_control(control_labels, column=column)
        self.qc_cells()
        self.quality_control(min_cells=min_cells)
        self.select_hvgs(n_top_genes=n_top_genes, column=column)
        self.normalize(target_sum=target_sum)
        self.log_transform()
        self.subset_genes()
        self.Y = self.mean_expressions(column=column)
        return self

    def _require(self, column):
        if column not in self.adata.obs.columns:
            raise KeyError(f"Column '{column}' not found in adata.obs")


def group_means(X, labels, var_names):
    """Mean of rows of X per label, via a sparse indicator matrix (no dense copy of X)."""
    labels = pd.Categorical(np.asarray(labels).astype(str))
    groups = labels.categories
    n = len(labels)

    indicator = sp.csr_matrix(
        (np.ones(n, dtype=np.float64), (labels.codes, np.arange(n))),
        shape=(len(groups), n),
    )
    counts = np.asarray(indicator.sum(axis=1)).ravel()
    sums = indicator @ X
    sums = sums.toarray() if sp.issparse(sums) else np.asarray(sums)

    return pd.DataFrame(sums / counts[:, None], index=groups, columns=var_names)
