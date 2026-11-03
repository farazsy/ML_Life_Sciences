from preprocessing.preprocessing import PrePipeline

# Initialize pipeline
pipeline = PrePipeline("AdamsonWeissman2016_GSM2406681_10X010.h5ad")

# Run preprocessing
adata = pipeline.run_pipeline(
    min_cells=5, min_counts=200, normalization="target", target_sum=1e6
)

# Display results
print("\nPreprocessing completed successfully!")

print(f"Number of cells: {adata.n_obs}")
print(f"Number of genes: {adata.n_vars}")

print("\nPerturbation groups:")
print(adata.obs["perturbation"].value_counts())

print("\nFirst five cells:")
print(adata.obs.head())

print("\nExpression matrix shape:")
print(adata.X.shape)

# Save processed dataset
adata.write_h5ad("processed_dataset.h5ad")

print("\nProcessed dataset saved successfully!")
