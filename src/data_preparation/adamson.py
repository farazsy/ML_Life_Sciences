"""
Adamson et al. (2016) specific preprocessing helpers (scPerturb raw release,
AdamsonWeissman2016_GSM2406681_10X010.h5ad).

Labels in obs['perturbation'] are guide-level, e.g. 'OST4_pDS353'. Controls are
'63(mod)_pBA580' and 'Gal4-4(mod)_pBA582' (pooled) and '62(mod)_pBA581' (2 cells, dropped).
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp

from build_target import group_means

CTRL = "ctrl"
CONTROL_GUIDES = ("63(mod)_pBA580", "Gal4-4(mod)_pBA582")
DROP_GUIDES = ("62(mod)_pBA581",)

# Old (hg19 / 10x 1.2) symbols of perturbed genes -> current HGNC symbols.
GENE_ALIASES = {
    "AARS": "AARS1",
    "ATP5B": "ATP5F1B",
    "C7orf26": "INTS15",
    "CARS": "CARS1",
    "DARS": "DARS1",
    "HARS": "HARS1",
    "MARS": "MARS1",
    "QARS": "QARS1",
    "SARS": "SARS1",
    "SLMO2": "PRELID3B",
    "TARS": "TARS1",
}


# Symbols that scPerturb de-duplicated onto the wrong Ensembl ID: symbol -> correct Ensembl ID.
# 'IER3IP1' in var is ENSG00000267228 (0 counts); the real IER3IP1 is stored as 'IER3IP1-1'.
SYMBOL_BY_ENSEMBL = {"IER3IP1": "ENSG00000134049"}


def fix_var_symbols(adata, fixes=SYMBOL_BY_ENSEMBL):
    """Give each symbol in `fixes` to its correct Ensembl ID; rename the wrong row to its Ensembl ID."""
    names = adata.var_names.to_numpy().astype(object)
    ens = adata.var["ensembl_id"].to_numpy()
    for symbol, ensembl in fixes.items():
        names[(names == symbol) & (ens != ensembl)] = None
        wrong = np.flatnonzero(names == None)  # noqa: E711
        names[wrong] = ens[wrong]
        names[ens == ensembl] = symbol
    adata.var_names = pd.Index(names.astype(str))
    if not adata.var_names.is_unique:
        raise ValueError("var_names not unique after symbol fixes")
    return adata


def parse_guide_gene(label):
    """'OST4_pDS353' -> 'OST4'; control guides -> 'ctrl'."""
    if label in CONTROL_GUIDES:
        return CTRL
    return label.rsplit("_", 1)[0]


def annotate_labels(adata, column="perturbation"):
    """
    Add obs['guide'] and obs['gene'] and drop the unused control guide.
    Expects NaN / '*' cells to be removed already (BuildTarget.filter_cells).
    """
    adata = adata[~adata.obs[column].isin(DROP_GUIDES).values].copy()
    guides = adata.obs[column].astype(str)
    adata.obs["guide"] = pd.Categorical(guides)
    adata.obs["gene"] = pd.Categorical(guides.map(parse_guide_gene))
    print(
        f"Labels: {adata.obs['guide'].nunique()} guides, "
        f"{adata.obs['gene'].nunique() - 1} genes + ctrl"
    )
    return adata


def _log1p_sparse(X):
    X = sp.csr_matrix(X, copy=True)
    X.data = np.log1p(X.data)
    return X


def guide_qc(
    adata,
    max_kd_ratio=0.7,
    min_agreement_r=0.5,
    min_ctrl_log_mean=0.2,
):
    """
    Per-guide quality table computed on CP10k data from raw counts.

    kd_ratio     mean target expression in guide cells / mean in control cells
                 (CRISPRi should knock it down; > max_kd_ratio => failed knockdown)
    agreement_r  for genes with >=2 guides: Pearson r between this guide's
                 delta-from-control profile (log1p CP10k pseudo-bulk) and the pooled
                 delta of the gene's other guides, over well-expressed genes
                 excluding the target itself (< min_agreement_r => disagrees)

    Also returns the control-guide agreement r (63(mod) vs Gal4-4(mod)).
    """
    guide = adata.obs["guide"].astype(str).to_numpy()
    gene = adata.obs["gene"].astype(str).to_numpy()
    ctrl_mask = gene == CTRL

    # Library sizes over ALL genes; then work on a column subset to save memory:
    # genes well expressed in controls + the perturbed target genes.
    lib = np.asarray(adata.X.sum(axis=1)).ravel()
    scale = sp.diags(1e4 / lib)

    ctrl_norm = sp.diags(1e4 / lib[ctrl_mask]) @ sp.csr_matrix(adata.X[ctrl_mask])
    ctrl_log_all = np.asarray(_log1p_sparse(ctrl_norm).mean(axis=0)).ravel()
    targets = set(gene) - {CTRL}
    cols = (ctrl_log_all > min_ctrl_log_mean) | adata.var_names.isin(list(targets))
    var = adata.var_names[cols]

    norm = scale @ sp.csr_matrix(adata.X[:, cols])  # CP10k, column subset
    log_norm = _log1p_sparse(norm)

    # Pseudo-bulk per guide (controls kept per guide for the control check).
    lin_means = group_means(norm, guide, var)
    log_means = group_means(log_norm, guide, var)

    ctrl_lin = np.asarray(norm[ctrl_mask].mean(axis=0)).ravel()
    ctrl_log = np.asarray(log_norm[ctrl_mask].mean(axis=0)).ravel()

    expressed = ctrl_log > min_ctrl_log_mean
    var_index = pd.Index(var)

    # Control guides should look alike.
    c1, c2 = (log_means.loc[g].to_numpy()[expressed] for g in CONTROL_GUIDES)
    ctrl_r = float(np.corrcoef(c1, c2)[0, 1])

    guide_gene = (
        pd.DataFrame({"guide": guide, "gene": gene})
        .drop_duplicates()
        .set_index("guide")["gene"]
    )
    n_cells = pd.Series(guide).value_counts()

    rows = []
    for g, target in guide_gene.items():
        if target == CTRL:
            continue

        # Knockdown ratio of the target gene.
        if target in var_index:
            j = var_index.get_loc(target)
            kd = lin_means.loc[g].iloc[j] / ctrl_lin[j] if ctrl_lin[j] > 0 else np.nan
        else:
            j, kd = None, np.nan

        # Agreement with sibling guides of the same gene.
        siblings = [s for s, t in guide_gene.items() if t == target and s != g]
        r = np.nan
        if siblings:
            cols = expressed.copy()
            if j is not None:
                cols[j] = False
            own = log_means.loc[g].to_numpy()[cols] - ctrl_log[cols]
            # Pool siblings weighted by cell count.
            w = n_cells[siblings].to_numpy(dtype=float)
            other = (log_means.loc[siblings].to_numpy()[:, cols] * w[:, None]).sum(0) / w.sum()
            other = other - ctrl_log[cols]
            r = float(np.corrcoef(own, other)[0, 1])

        rows.append(
            {
                "guide": g,
                "gene": target,
                "n_cells": int(n_cells[g]),
                "n_guides_for_gene": len(siblings) + 1,
                "kd_ratio": kd,
                "agreement_r": r,
            }
        )

    qc = pd.DataFrame(rows).set_index("guide")
    qc["failed_kd"] = qc["kd_ratio"] > max_kd_ratio
    qc["disagrees"] = qc["agreement_r"] < min_agreement_r
    # Disagreement alone can't tell which of two guides is wrong (and weak-effect genes
    # give noisy deltas), so it only removes a guide whose knockdown can't be measured.
    qc["flagged"] = qc["failed_kd"] | (qc["disagrees"] & qc["kd_ratio"].isna())
    return qc.sort_values(["gene", "n_cells"], ascending=[True, False]), ctrl_r


def drop_guides(adata, guides):
    """Remove all cells carrying any of the given guides."""
    keep = ~adata.obs["guide"].isin(list(guides)).values
    print(f"Dropping {len(guides)} flagged guides ({(~keep).sum()} cells)")
    out = adata[keep].copy()
    out.obs["guide"] = out.obs["guide"].cat.remove_unused_categories()
    out.obs["gene"] = out.obs["gene"].cat.remove_unused_categories()
    return out


def filter_min_cells(adata, min_cells=50):
    """Keep perturbed genes with >= min_cells cells after merging guides (ctrl always kept)."""
    counts = adata.obs["gene"].value_counts()
    small = [g for g, n in counts.items() if n < min_cells and g != CTRL]
    keep = ~adata.obs["gene"].isin(small).values
    print(f"Dropping {len(small)} genes with < {min_cells} cells: {sorted(small)}")
    out = adata[keep].copy()
    out.obs["gene"] = out.obs["gene"].cat.remove_unused_categories()
    out.obs["guide"] = out.obs["guide"].cat.remove_unused_categories()
    return out, small


def harmonize_symbols(adata, aliases=GENE_ALIASES):
    """Rename old symbols in var_names and obs['gene']; keep originals in var['symbol_orig']."""
    clashes = [new for old, new in aliases.items() if old in adata.var_names and new in adata.var_names]
    if clashes:
        raise ValueError(f"Alias targets already present in var_names: {clashes}")

    adata.var["symbol_orig"] = adata.var_names.astype(str)
    adata.var_names = [aliases.get(g, g) for g in adata.var_names]
    adata.obs["gene"] = adata.obs["gene"].cat.rename_categories(
        lambda g: aliases.get(g, g)
    )
    renamed = sorted(set(aliases) & set(adata.var["symbol_orig"]))
    print(f"Renamed {len(renamed)} genes to current HGNC symbols: {renamed}")
    return adata
