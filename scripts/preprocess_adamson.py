"""
Preprocess the raw scPerturb Adamson (2016) file into a GenePert-ready AnnData.

    python scripts/preprocess_adamson.py \
        --input AdamsonWeissman2016_GSM2406681_10X010.h5ad \
        --output data/processed/adamson_genepert.h5ad

Steps (see GenePT_GenePert_study_notes.md, sections 5-7):
 1. load raw counts
 2. drop NaN / '*' labels, drop 62(mod), pool 63(mod) + Gal4-4(mod) as ctrl
 3. cell QC (mito > 15 %, < 500 genes, counts > median + 5 MAD); genes in < 10 cells
 4. guide QC: knockdown ratio and sibling-guide agreement -> drop flagged guides
 5. merge guides into genes; keep genes with >= 50 cells
 6. harmonize old gene symbols to current HGNC
 7. HVGs (seurat_v3, control cells only) + perturbed genes
 8. normalize to 1e4, log1p, subset genes (raw counts kept in layers['counts'])
 9. obs['condition'] = 'GENE+ctrl' / 'ctrl'
10. validate and write h5ad + guide_qc.csv + pseudo-bulk means
"""

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "data_preparation"))
sys.path.insert(0, str(ROOT / "src" / "data_ingestion"))

import adamson  # noqa: E402
from build_target import BuildTarget  # noqa: E402
from utils import read_target  # noqa: E402


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", default=str(ROOT / "AdamsonWeissman2016_GSM2406681_10X010.h5ad"))
    p.add_argument("--output", default=str(ROOT / "data" / "processed" / "adamson_genepert.h5ad"))
    p.add_argument("--n-hvg", type=int, default=5000)
    p.add_argument("--min-cells-per-pert", type=int, default=50)
    p.add_argument("--min-cells-per-gene", type=int, default=10, help="gene filter: detected in >= N cells")
    p.add_argument("--max-pct-mito", type=float, default=15.0)
    p.add_argument("--min-genes", type=int, default=500)
    p.add_argument("--max-counts-mad", type=float, default=5.0)
    p.add_argument("--max-kd-ratio", type=float, default=0.7)
    p.add_argument("--min-agreement-r", type=float, default=0.5)
    p.add_argument("--keep-failed-kd", action="store_true", help="keep flagged guides (report only)")
    p.add_argument("--no-harmonize", action="store_true", help="keep original (hg19) gene symbols")
    return p.parse_args(argv)


def validate(adata, min_cells_per_pert):
    cond = adata.obs["condition"].astype(str)
    assert cond.eq("ctrl").sum() >= 6000, "expected >= 6000 control cells"
    assert cond[cond != "ctrl"].str.endswith("+ctrl").all(), "bad condition labels"

    sizes = cond.value_counts()
    assert (sizes >= min_cells_per_pert).all(), "condition below min cells"

    genes = [c.replace("+ctrl", "") for c in sizes.index if c != "ctrl"]
    missing = [g for g in genes if g not in adata.var_names]
    assert not missing, f"perturbed genes missing from var: {missing}"

    X = adata.X
    data = X.data if hasattr(X, "data") else np.asarray(X)
    assert np.isfinite(data).all(), "non-finite values in X"
    assert data.max() < 15, "X does not look log-transformed"
    return genes


def main(argv=None):
    args = parse_args(argv)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # 1. load
    adata = read_target(args.input)
    sample = adata.X[:1000].toarray() if hasattr(adata.X, "toarray") else adata.X[:1000]
    assert np.allclose(sample, np.round(sample)), "expected raw integer counts in X"
    print(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")
    adata = adamson.fix_var_symbols(adata)

    # 2. labels
    bt = BuildTarget(adata)
    bt.filter_cells(column="perturbation", drop_labels=("*",))
    bt.adata = adamson.annotate_labels(bt.adata)

    # 3. QC
    bt.qc_cells(max_pct_mito=args.max_pct_mito, min_genes=args.min_genes, max_counts_mad=args.max_counts_mad)
    targets = set(bt.adata.obs["gene"].astype(str)) - {adamson.CTRL}
    bt.quality_control(min_cells=args.min_cells_per_gene, protect_genes=targets)

    # 4. guide QC
    qc, ctrl_r = adamson.guide_qc(
        bt.adata, max_kd_ratio=args.max_kd_ratio, min_agreement_r=args.min_agreement_r
    )
    print(f"Control guides agreement r = {ctrl_r:.4f}")
    if ctrl_r < 0.98:
        print("Warning: control guides differ more than expected")

    flagged = qc.index[qc["flagged"]].tolist()
    print(f"Flagged guides ({len(flagged)}):")
    print(qc.loc[flagged, ["gene", "n_cells", "kd_ratio", "agreement_r"]].round(3).to_string())
    qc["kept"] = ~qc["flagged"] | args.keep_failed_kd
    if not args.keep_failed_kd and flagged:
        bt.adata = adamson.drop_guides(bt.adata, flagged)

    # 5. merge guides -> genes, min cells
    bt.adata, small = adamson.filter_min_cells(bt.adata, min_cells=args.min_cells_per_pert)
    qc.loc[qc["gene"].isin(small), "kept"] = False

    # 6. symbols
    if not args.no_harmonize:
        bt.adata = adamson.harmonize_symbols(bt.adata)
        qc["gene_hgnc"] = qc["gene"].map(lambda g: adamson.GENE_ALIASES.get(g, g))

    # 7-8. HVGs on raw counts, then normalize over all genes, then subset
    perturbed = [g for g in bt.adata.obs["gene"].cat.categories if g != adamson.CTRL]
    bt.select_hvgs(n_top_genes=args.n_hvg, column="gene", ctrl_label=adamson.CTRL, add_genes=perturbed)
    bt.normalize(target_sum=1e4)
    bt.log_transform()
    bt.subset_genes(keep_counts_layer=True)

    # 9. condition labels
    bt.set_condition(gene_column="gene", ctrl_label=adamson.CTRL)
    bt.adata.uns["preprocessing"] = {k: (str(v) if v is None else v) for k, v in vars(args).items()}

    # 10. validate + save
    genes = validate(bt.adata, args.min_cells_per_pert)
    means = bt.mean_expressions(column="condition")

    bt.adata.write_h5ad(out_path)
    qc_path = out_path.with_name("guide_qc.csv")
    means_path = out_path.with_name(out_path.stem + "_pseudobulk.csv")
    qc.to_csv(qc_path)
    means.to_csv(means_path)

    a = bt.adata
    print("\n=== Summary ===")
    print(f"Cells: {a.n_obs} ({(a.obs['condition'] == 'ctrl').sum()} ctrl)")
    print(f"Genes: {a.n_vars}")
    print(f"Perturbations: {len(genes)}")
    print(a.obs["condition"].value_counts().describe().round(1).to_string())
    print(f"Wrote {out_path}\n      {qc_path}\n      {means_path}")
    return a


if __name__ == "__main__":
    main()
