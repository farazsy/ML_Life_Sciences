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

CTRL = "ctrl"


class BuildTarget:
    """Turn a single-cell perturbation AnnData into GenePert-ready cells and pseudo-bulk targets."""

    def __init__(self, adata=None, copy=False):
        """Create the builder, optionally loading adata now (in place unless copy=True)."""
        # Copying a large raw matrix doubles memory, so default to in-place.
        self.adata = None
        self.Y = None
        self.targets = None
        self.keep_genes = None
        self.is_counts = None
        self._counts = None
        if adata is not None:
            self.set_data(adata, copy=copy)

    def set_data(self, adata, copy=False):
        """Attach an AnnData, make gene names unique and detect whether X is raw counts."""
        self.adata = adata.copy() if copy else adata
        if not self.adata.var_names.is_unique:
            n_dup = int(self.adata.var_names.duplicated().sum())
            print(f"Warning: {n_dup} duplicated gene names; making them unique (GENE-1, ...)")
            self.adata.var_names_make_unique()
        self.is_counts = is_raw_counts(self.adata.X)
        print(
            f"Input: {self.adata.n_obs} cells x {self.adata.n_vars} genes, "
            f"X looks like {'raw counts' if self.is_counts else 'normalized values'}"
        )
        return self

    # ------------------------------------------------------------------ labels
    def filter_cells(self, column="perturbation", drop_labels=("*",)):
        """Remove cells without a valid perturbation label (NaN, empty or in drop_labels)."""
        self._require(column)

        labels = self.adata.obs[column]
        mask = labels.notna() & ~labels.isin(list(drop_labels))
        mask &= labels.astype(str).str.strip() != ""

        self.adata = self.adata[mask.values].copy()
        self.adata.obs[column] = self.adata.obs[column].astype(str)

        print(f"Cells after label filtering: {self.adata.n_obs}")
        return self

    def format_control(self, control_labels, column="perturbation", ctrl_label=CTRL):
        """Map all control labels (non-targeting guides, 'control', ...) to one control label."""
        self._require(column)

        labels = self.adata.obs[column].astype(str)
        control_labels = list(control_labels)
        present = [c for c in control_labels if c in set(labels)]
        missing = [c for c in control_labels if c not in set(labels)]
        if missing:
            print(f"Warning: control labels not found in obs['{column}']: {missing}")
        if not present and ctrl_label not in set(labels):
            raise ValueError(
                f"No control cells: none of {control_labels} (or '{ctrl_label}') in obs['{column}']"
            )

        self.adata.obs[column] = labels.replace(present, ctrl_label).astype("category")
        print(f"Control cells: {(self.adata.obs[column] == ctrl_label).sum()}")
        return self

    def annotate_genes(
        self,
        column="perturbation",
        label_to_gene=None,
        combo_sep="+",
        ctrl_label=CTRL,
        out_column="gene",
    ):
        """
        Derive obs[out_column] = perturbed gene(s) in canonical form ('ctrl', 'G', 'G1+G2').

        label_to_gene: optional callable for guide-level labels, e.g.
            lambda s: s.rsplit("_", 1)[0]   # 'OST4_pDS353' -> 'OST4'
        It is applied to every non-control label before splitting combos on combo_sep.
        """
        self._require(column)

        labels = self.adata.obs[column].astype(str)
        mapping = {}
        for label in labels.unique():
            if label == ctrl_label:
                mapping[label] = ctrl_label
                continue
            gene = label_to_gene(label) if label_to_gene is not None else label
            mapping[label] = canonical_gene_label(gene, ctrl_label=ctrl_label, combo_sep=combo_sep)

        self.adata.obs[out_column] = pd.Categorical(labels.map(mapping))
        n_pert = self.adata.obs[out_column].nunique() - int(ctrl_label in mapping.values())
        n_combo = sum("+" in g for g in set(mapping.values()))
        print(f"Perturbations: {n_pert} ({n_combo} combinatorial) + {ctrl_label}")
        return self

    def perturbed_genes(self, gene_column="gene", ctrl_label=CTRL):
        """All individual genes targeted by any perturbation (combos are split)."""
        self._require(gene_column)
        genes = set()
        for label in self.adata.obs[gene_column].astype(str).unique():
            genes.update(g for g in label.split("+") if g and g != ctrl_label)
        return sorted(genes)

    def filter_min_cells(self, min_cells=50, gene_column="gene", ctrl_label=CTRL):
        """Drop perturbations with fewer than min_cells cells (control is always kept)."""
        self._require(gene_column)

        sizes = self.adata.obs[gene_column].value_counts()
        small = sorted(g for g, n in sizes.items() if n < min_cells and g != ctrl_label)
        if small:
            print(f"Dropping {len(small)} perturbations with < {min_cells} cells: {small}")
            keep = ~self.adata.obs[gene_column].isin(small).values
            self.adata = self.adata[keep].copy()
            for col in self.adata.obs.select_dtypes("category"):
                self.adata.obs[col] = self.adata.obs[col].cat.remove_unused_categories()
        return self

    # ---------------------------------------------------------------------- QC
    def qc_cells(self, max_pct_mito=15.0, min_genes=500, max_counts_mad=5.0, mito_prefix="MT-"):
        """
        Cell-level QC on raw counts. Set a threshold to None to skip that filter.

        - drop cells with percent_mito > max_pct_mito (genes starting with mito_prefix,
          case-insensitive, so human 'MT-' and mouse 'mt-' both match)
        - drop cells with fewer than min_genes detected genes
        - drop likely doublets: total counts > median + max_counts_mad * MAD
        """
        self._require_counts("qc_cells")
        X = self.adata.X
        n_counts = np.asarray(X.sum(axis=1)).ravel()
        n_genes = np.asarray((X > 0).sum(axis=1)).ravel()
        mask = n_counts > 0
        report = [f"empty: {(~mask).sum()}"]

        if max_pct_mito is not None:
            if "percent_mito" in self.adata.obs:
                pct_mito = self.adata.obs["percent_mito"].to_numpy(dtype=float)
            else:
                mito = self.adata.var_names.str.upper().str.startswith(mito_prefix.upper())
                if not mito.any():
                    print(f"Note: no genes start with '{mito_prefix}'; skipping mito filter")
                mito_counts = np.asarray(X[:, mito].sum(axis=1)).ravel()
                pct_mito = 100 * np.divide(
                    mito_counts, n_counts, out=np.zeros_like(n_counts, dtype=float), where=n_counts > 0
                )
            mask &= pct_mito <= max_pct_mito
            report.append(f"mito>{max_pct_mito}%: {(pct_mito > max_pct_mito).sum()}")

        if min_genes is not None:
            mask &= n_genes >= min_genes
            report.append(f"genes<{min_genes}: {(n_genes < min_genes).sum()}")

        if max_counts_mad is not None:
            median = np.median(n_counts)
            mad = np.median(np.abs(n_counts - median))
            max_counts = median + max_counts_mad * mad
            mask &= n_counts <= max_counts
            report.append(f"counts>{max_counts:.0f}: {(n_counts > max_counts).sum()}")

        print(f"Cell QC: {', '.join(report)} -> removing {(~mask).sum()} cells")
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
            self._require_counts("quality_control(min_counts=...)")
            sc.pp.filter_cells(self.adata, min_counts=min_counts)

        print(f"After QC: {self.adata.n_obs} cells, {self.adata.n_vars} genes")
        return self

    # ----------------------------------------------------------- HVG selection
    def select_hvgs(
        self,
        n_top_genes=5000,
        column="gene",
        ctrl_label=CTRL,
        control_only=True,
        add_genes=(),
        flavor=None,
    ):
        """
        Choose HVGs and remember them (plus ``add_genes``) for ``subset_genes``.

        With control_only=True HVGs are computed on control cells only, so the gene set
        does not depend on which perturbations end up in train or test folds.
        flavor=None picks 'seurat_v3' for raw counts (must run before ``normalize``)
        and 'seurat' for log-normalized data. The subset itself is applied later by
        ``subset_genes`` so that normalization uses total counts over all genes.
        """
        self._require(column)
        if flavor is None:
            flavor = "seurat_v3" if self.is_counts else "seurat"
        if flavor == "seurat_v3":
            self._require_counts("select_hvgs(flavor='seurat_v3')")

        if control_only:
            is_ctrl = (self.adata.obs[column] == ctrl_label).values
            if is_ctrl.sum() < 2:
                raise ValueError(f"control_only=True but obs['{column}'] has no '{ctrl_label}' cells")
            ref = self.adata[is_ctrl].copy()
        else:
            ref = self.adata.copy()

        n_top = min(n_top_genes, ref.n_vars)
        sc.pp.highly_variable_genes(ref, n_top_genes=n_top, flavor=flavor)
        hvgs = ref.var_names[ref.var["highly_variable"]]

        hvg_set = set(hvgs)
        extra = [g for g in add_genes if g in self.adata.var_names and g not in hvg_set]
        missing = [g for g in add_genes if g not in self.adata.var_names]
        if missing:
            print(f"Warning: {len(missing)} add_genes not in var: {missing}")

        keep = hvg_set | set(extra)
        self.keep_genes = [g for g in self.adata.var_names if g in keep]  # keep var order
        self.adata.var["highly_variable"] = self.adata.var_names.isin(hvgs)

        # Stash raw counts of the kept genes for layers['counts'] (only if X is counts).
        self._counts = self.adata[:, self.keep_genes].X.copy() if self.is_counts else None

        print(f"HVGs ({flavor}): {len(hvgs)} + {len(extra)} added genes = {len(self.keep_genes)}")
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
        self._require_counts("normalize")
        sc.pp.normalize_total(self.adata, target_sum=target_sum)
        self.is_counts = False
        return self

    def log_transform(self):
        """Apply log1p transformation."""
        if self.is_counts:
            raise RuntimeError("log_transform on raw counts: call normalize first.")
        sc.pp.log1p(self.adata)
        return self

    # ---------------------------------------------------------- conditions
    def set_condition(self, gene_column="gene", ctrl_label=CTRL, out_column="condition", combo_sep="+"):
        """Build GEARS/GenePert labels: 'ctrl', 'GENE+ctrl', 'GENE1+GENE2'."""
        self._require(gene_column)

        genes = self.adata.obs[gene_column].astype(str)
        mapping = {
            g: condition_label(g, ctrl_label=ctrl_label, combo_sep=combo_sep) for g in genes.unique()
        }
        self.adata.obs[out_column] = pd.Categorical(genes.map(mapping))
        return self

    def mean_expressions(self, column="condition"):
        """Pseudo-bulk: mean expression per condition (conditions x genes DataFrame)."""
        self._require(column)
        return group_means(self.adata.X, self.adata.obs[column], self.adata.var_names)

    def build_targets(self, column="condition", ctrl_label=CTRL):
        """
        Training targets: one row per perturbation (control excluded).

        Returns an AnnData (perturbations x genes) with
            obs_names        perturbed gene(s): 'GENE' or 'GENE1+GENE2'
            obs['condition'] original condition label ('GENE+ctrl', ...)
            obs['n_cells']   cells averaged into the row
            X                mean expression per perturbation
            layers['delta']  X - control mean
            var['ctrl_mean'] mean expression of control cells
        """
        self._require(column)

        means = self.mean_expressions(column=column)
        if ctrl_label not in means.index:
            raise ValueError(f"No '{ctrl_label}' rows in obs['{column}']")

        ctrl_mean = means.loc[ctrl_label].to_numpy()
        pert = means.drop(index=ctrl_label)
        sizes = self.adata.obs[column].astype(str).value_counts()

        names = [gene_from_condition(c, ctrl_label=ctrl_label) for c in pert.index]
        X = pert.to_numpy()
        targets = ad.AnnData(
            X=X,
            obs=pd.DataFrame(
                {"condition": pert.index.astype(str), "n_cells": sizes[pert.index].to_numpy()},
                index=pd.Index(names),
            ),
            var=pd.DataFrame({"ctrl_mean": ctrl_mean}, index=pert.columns.copy()),
            layers={"delta": X - ctrl_mean},
        )
        if not targets.obs_names.is_unique:
            raise ValueError("Duplicate perturbation names in targets (e.g. 'A+B' and 'B+A')")

        self.Y = means
        self.targets = targets
        print(f"Targets: {targets.n_obs} perturbations x {targets.n_vars} genes")
        return targets

    # ---------------------------------------------------------- validation
    def validate(self, column="condition", ctrl_label=CTRL, min_cells=None, require_targets_in_var=False):
        """
        Check the cell-level output contract; raises ValueError on failure.
        Returns the list of perturbation names (conditions without the '+ctrl' suffix).
        """
        self._require(column)
        errors = []
        cond = self.adata.obs[column].astype(str)

        if not (cond == ctrl_label).any():
            errors.append(f"no '{ctrl_label}' cells")
        bad = sorted(c for c in cond.unique() if c != ctrl_label and "+" not in c)
        if bad:
            errors.append(f"labels not in 'GENE+{ctrl_label}' / 'G1+G2' form: {bad[:5]}")

        sizes = cond.value_counts()
        if min_cells is not None:
            small = sorted(sizes.index[(sizes < min_cells) & (sizes.index != ctrl_label)])
            if small:
                errors.append(f"conditions with < {min_cells} cells: {small[:5]}")

        perts = [gene_from_condition(c, ctrl_label=ctrl_label) for c in sizes.index if c != ctrl_label]
        if require_targets_in_var:
            genes = {g for p in perts for g in p.split("+")}
            missing = sorted(genes - set(self.adata.var_names))
            if missing:
                errors.append(f"perturbed genes missing from var: {missing[:5]}")

        X = self.adata.X
        data = X.data if sp.issparse(X) else np.asarray(X)
        if not np.isfinite(data).all():
            errors.append("non-finite values in X")
        elif data.size and data.max() >= 15:
            errors.append(f"X max {data.max():.1f} does not look log-transformed")

        if errors:
            raise ValueError("Validation failed:\n  - " + "\n  - ".join(errors))
        print(f"Validation passed: {len(perts)} perturbations, {sizes.get(ctrl_label, 0)} control cells")
        return perts

    # ---------------------------------------------------------- pipeline
    def run_pipeline(
        self,
        adata=None,
        control_labels=(),
        column="perturbation",
        label_to_gene=None,
        combo_sep="+",
        drop_labels=("*",),
        ctrl_label=CTRL,
        max_pct_mito=15.0,
        min_genes=500,
        max_counts_mad=5.0,
        min_cells=10,
        min_cells_per_pert=50,
        n_top_genes=5000,
        hvg_control_only=True,
        target_sum=1e4,
    ):
        """
        Generic pipeline: raw (or log-normalized) AnnData -> GenePert-ready AnnData + targets.

        Dataset-specific choices are arguments:
            control_labels   labels of control cells, e.g. ('non-targeting',)
            label_to_gene    guide label -> gene symbol, e.g. lambda s: s.rsplit('_', 1)[0]
            combo_sep        separator of combinatorial labels in the input ('+', '_', ...)
        Set a QC threshold / n_top_genes / min_cells_per_pert to None to skip that step.
        If X is not raw counts, cell QC and normalization are skipped.
        """
        if adata is not None:
            self.set_data(adata)
        if self.adata is None:
            raise ValueError("No data: pass adata to BuildTarget(...) or run_pipeline(...)")

        self.filter_cells(column=column, drop_labels=drop_labels)
        self.format_control(control_labels, column=column, ctrl_label=ctrl_label)
        self.annotate_genes(column=column, label_to_gene=label_to_gene, combo_sep=combo_sep, ctrl_label=ctrl_label)
        targets = self.perturbed_genes(ctrl_label=ctrl_label)

        if self.is_counts:
            self.qc_cells(max_pct_mito=max_pct_mito, min_genes=min_genes, max_counts_mad=max_counts_mad)
        else:
            print("X is not raw counts: skipping cell QC and normalization (assuming log-normalized)")
        self.quality_control(min_cells=min_cells, protect_genes=targets)
        if min_cells_per_pert:
            self.filter_min_cells(min_cells=min_cells_per_pert, ctrl_label=ctrl_label)
            targets = self.perturbed_genes(ctrl_label=ctrl_label)

        if n_top_genes:
            self.select_hvgs(
                n_top_genes=n_top_genes,
                column="gene",
                ctrl_label=ctrl_label,
                control_only=hvg_control_only,
                add_genes=targets,
            )
        if self.is_counts:
            self.normalize(target_sum=target_sum)
            self.log_transform()
        if n_top_genes:
            self.subset_genes()

        self.set_condition(gene_column="gene", ctrl_label=ctrl_label)
        self.validate(ctrl_label=ctrl_label, min_cells=min_cells_per_pert or None)
        self.build_targets(ctrl_label=ctrl_label)
        return self

    # ---------------------------------------------------------- helpers
    def _require(self, column):
        """Raise if no data is loaded or obs[column] does not exist."""
        if self.adata is None:
            raise ValueError("No data loaded")
        if column not in self.adata.obs.columns:
            raise KeyError(f"Column '{column}' not found in adata.obs")

    def _require_counts(self, step):
        """Raise if `step` needs raw counts but X is already normalized."""
        if not self.is_counts:
            raise RuntimeError(
                f"{step} needs raw counts, but X looks normalized/log-transformed. "
                "Skip this step for pre-normalized data."
            )


def is_raw_counts(X, n_rows=1000):
    """True if a sample of X is non-negative and integer-valued."""
    sample = X[: min(n_rows, X.shape[0])]
    values = sample.data if sp.issparse(sample) else np.asarray(sample).ravel()
    if values.size == 0:
        return True
    return bool((values >= 0).all() and np.allclose(values, np.round(values)))


def canonical_gene_label(label, ctrl_label=CTRL, combo_sep="+"):
    """
    Normalize a perturbation label to gene form: drop control parts, join combos with '+'.

    'G1' -> 'G1';  'G1+ctrl' / 'ctrl+G1' -> 'G1';  'G1_G2' (combo_sep='_') -> 'G1+G2';
    'ctrl' / 'ctrl+ctrl' -> 'ctrl'.
    """
    parts = [p.strip() for p in str(label).split(combo_sep)]
    genes = [p for p in parts if p and p != ctrl_label]
    return "+".join(genes) if genes else ctrl_label


def condition_label(gene_label, ctrl_label=CTRL, combo_sep="+"):
    """Gene label -> condition: 'ctrl' -> 'ctrl', 'G' -> 'G+ctrl', 'G1+G2' -> 'G1+G2'."""
    gene = canonical_gene_label(gene_label, ctrl_label=ctrl_label, combo_sep=combo_sep)
    if gene == ctrl_label or "+" in gene:
        return gene
    return f"{gene}+{ctrl_label}"


def gene_from_condition(condition, ctrl_label=CTRL):
    """Condition -> perturbation name: 'G+ctrl' -> 'G', 'G1+G2' -> 'G1+G2'."""
    return canonical_gene_label(condition, ctrl_label=ctrl_label, combo_sep="+")


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
