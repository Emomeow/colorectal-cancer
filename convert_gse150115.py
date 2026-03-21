#!/usr/bin/env python3
"""
Convert GSE150115 to h5ad format
Data structure: genes × cells (rows=genes, columns=cells)
"""

import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix, vstack
import anndata as ad
import glob
import os
import gc
from datetime import datetime

print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Starting GSE150115 conversion")
print("=" * 80)

# Step 1: Get all sample files
print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Step 1: Finding sample files...")
sample_files = sorted(glob.glob('data/GSE150115_extracted/*.umi.dge*.txt.gz'))
n_samples = len(sample_files)
print(f"Found {n_samples} DGE sample files")

if n_samples == 0:
    raise ValueError("No sample files found!")

for i, f in enumerate(sample_files, 1):
    print(f"  {i}. {os.path.basename(f)}")

# Step 2: Process each sample
print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Step 2: Processing samples...")
all_sparse_matrices = []
all_cell_ids = []
all_genes = None
total_cells = 0

for i, sample_file in enumerate(sample_files, 1):
    sample_name = os.path.basename(sample_file).split('_')[0]  # Extract GSM ID
    print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Processing {i}/{n_samples}: {sample_name}")

    # Read sample file (genes × cells)
    # First column is GENE, rest are cell barcodes
    sample_df = pd.read_csv(sample_file, compression='gzip', sep='\t', index_col=0)

    n_genes_in_sample = sample_df.shape[0]
    n_cells_in_sample = sample_df.shape[1]
    print(f"  Sample dimensions: {n_genes_in_sample} genes × {n_cells_in_sample} cells")

    # Store gene names from first sample
    if all_genes is None:
        all_genes = sample_df.index.tolist()
        print(f"  Using gene list from first sample: {len(all_genes)} genes")
    else:
        # Verify all samples have same genes
        sample_genes = sample_df.index.tolist()
        if sample_genes != all_genes:
            print(f"  WARNING: Gene list differs from first sample!")
            # Align to first sample's gene list
            sample_df = sample_df.reindex(all_genes, fill_value=0)
            print(f"  Aligned to reference gene list")

    # Transpose to cells × genes
    sample_df_t = sample_df.T

    # Convert to sparse matrix
    print(f"  Converting to sparse matrix...")
    sample_sparse = csr_matrix(sample_df_t.values, dtype=np.float32)

    sparsity = 100 * (1 - sample_sparse.nnz / (sample_sparse.shape[0] * sample_sparse.shape[1]))
    print(f"  Sparse matrix: {sample_sparse.shape[0]:,} cells × {sample_sparse.shape[1]:,} genes")
    print(f"  Non-zero: {sample_sparse.nnz:,}, Sparsity: {sparsity:.1f}%")

    all_sparse_matrices.append(sample_sparse)

    # Add sample prefix to cell IDs to avoid conflicts
    sample_cell_ids = [f"{sample_name}_{cell_id}" for cell_id in sample_df_t.index]
    all_cell_ids.extend(sample_cell_ids)

    total_cells += n_cells_in_sample

    # Clean up
    del sample_df, sample_df_t, sample_sparse
    gc.collect()

    print(f"  Accumulated cells: {total_cells:,}")

# Step 3: Vertically stack all samples
print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Step 3: Stacking all samples...")
print(f"Stacking {len(all_sparse_matrices)} sparse matrices...")
X_sparse = vstack(all_sparse_matrices, format='csr', dtype=np.float32)

overall_sparsity = 100 * (1 - X_sparse.nnz / (X_sparse.shape[0] * X_sparse.shape[1]))
print(f"Final matrix: {X_sparse.shape[0]:,} cells × {X_sparse.shape[1]:,} genes")
print(f"Non-zero elements: {X_sparse.nnz:,}")
print(f"Sparsity: {overall_sparsity:.2f}%")

# Clean up intermediate matrices
del all_sparse_matrices
gc.collect()

# Step 4: Create AnnData object
print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Step 4: Creating AnnData object...")

# Create obs (cell metadata)
obs_df = pd.DataFrame(index=all_cell_ids)
obs_df['dataset'] = 'GSE150115'

# Create var (gene metadata)
var_df = pd.DataFrame(index=all_genes)

# Create AnnData
adata = ad.AnnData(
    X=X_sparse,
    obs=obs_df,
    var=var_df
)

print(f"AnnData created: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
print(f"  X type: {type(adata.X)}")
print(f"  X dtype: {adata.X.dtype}")

# Verify data integrity
print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Step 5: Verifying data...")
print(f"  Total counts per cell (first 3):")
for i in range(min(3, adata.n_obs)):
    cell_counts = adata.X[i, :].sum()
    print(f"    Cell {i+1}: {cell_counts:,.0f} counts")

# Check for negative values (should be none for UMI counts)
min_val = adata.X.min()
max_val = adata.X.max()
print(f"  Value range: [{min_val:.0f}, {max_val:.0f}]")
if min_val < 0:
    print(f"  WARNING: Negative values detected! This may not be raw UMI counts.")
else:
    print(f"  OK: All values non-negative (raw UMI counts)")

# Step 6: Save to h5ad
print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Step 6: Saving to h5ad...")
output_dir = 'processed_data'
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, 'GSE150115.h5ad')

adata.write_h5ad(output_file, compression='gzip')

file_size_mb = os.path.getsize(output_file) / (1024**2)
print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] SUCCESS!")
print(f"Output: {output_file}")
print(f"File size: {file_size_mb:.1f} MB")
print(f"Final dimensions: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
print(f"Sparsity: {overall_sparsity:.2f}%")
print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Conversion complete!")
