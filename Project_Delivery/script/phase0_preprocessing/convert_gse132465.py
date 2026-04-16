"""
Convert GSE132465 to h5ad format
Larger dataset: 63,689 cells
Use SPARSE matrix from the start to avoid memory issues
"""

import pandas as pd
import scanpy as sc
import anndata as ad
from pathlib import Path
import sys
from scipy.sparse import csr_matrix
import numpy as np

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8')

data_dir = Path("data")
output_dir = Path("processed_data")
output_dir.mkdir(exist_ok=True)

print("="*80)
print("Converting GSE132465 to h5ad format")
print("="*80)

# Load annotation first (small file)
print("\nStep 1: Loading annotation...")
annot_file = data_dir / "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
annotation = pd.read_csv(annot_file, compression='gzip', sep='\t', index_col=0)
print(f"✓ Loaded {len(annotation)} cells")
print(f"  Columns: {list(annotation.columns)}")
print(f"  First cell ID: {annotation.index[0]}")

# Load count matrix WITH CHUNKS to avoid memory issues
print("\nStep 2: Loading count matrix (LARGE FILE - may take time)...")
count_file = data_dir / "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"
print(f"  File size: {count_file.stat().st_size / 1024 / 1024:.1f} MB")

# Strategy: Read in chunks if needed, but try full read first since we convert to sparse immediately
print(f"  Reading matrix (genes x cells)...")
count_matrix = pd.read_csv(count_file, compression='gzip', sep='\t', index_col=0)
print(f"✓ Loaded matrix: {count_matrix.shape[0]} genes x {count_matrix.shape[1]} cells")
print(f"  First cell ID in matrix: {count_matrix.columns[0]}")
print(f"  First gene: {count_matrix.index[0]}")

# Check data type
print(f"\nStep 3: Checking data type...")
print(f"  Matrix dtype: {count_matrix.dtypes[0]}")
print(f"  Sample values:\n{count_matrix.iloc[:3, :3]}")

# Match cells (no conversion needed for this dataset)
print(f"\nStep 4: Matching cells...")
common_cells = list(set(count_matrix.columns) & set(annotation.index))
print(f"  Common cells: {len(common_cells)} / {len(count_matrix.columns)}")

if len(common_cells) == len(count_matrix.columns):
    print(f"  ✅ Perfect match!")
elif len(common_cells) < len(count_matrix.columns):
    print(f"  ⚠ Warning: {len(count_matrix.columns) - len(common_cells)} cells not matched")

# Create AnnData object with SPARSE matrix
print(f"\nStep 5: Creating AnnData object with SPARSE matrix...")

# Transpose to cells x genes and select common cells
print(f"  Transposing and selecting matched cells...")
dense_matrix = count_matrix.T.loc[common_cells, :]
print(f"  Shape: {dense_matrix.shape[0]} cells x {dense_matrix.shape[1]} genes")

# Convert to sparse
print(f"  Converting to sparse matrix (this may take a moment)...")
sparse_matrix = csr_matrix(dense_matrix.values)

# Calculate sparsity
print(f"  Calculating sparsity...")
n_zeros = (dense_matrix.values == 0).sum()
total = dense_matrix.values.size
sparsity = (n_zeros / total) * 100
print(f"  Sparsity: {sparsity:.1f}% zeros")

# Clear the dense matrix to free memory
del dense_matrix
print(f"  ✓ Freed memory from dense matrix")

# Create AnnData
print(f"\nStep 6: Creating AnnData object...")
adata = ad.AnnData(X=sparse_matrix)
adata.obs_names = common_cells
adata.var_names = count_matrix.index
adata.obs = annotation.loc[common_cells, :]
adata.obs['dataset'] = 'GSE132465'
adata.obs['batch'] = 'GSE132465'

# Store raw counts (also sparse)
adata.layers['counts'] = adata.X.copy()

print(f"✓ Created AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes")
print(f"  X type: {type(adata.X)}")
print(f"  X dtype: {adata.X.dtype}")
print(f"  Obs columns: {list(adata.obs.columns)}")

# Verify data
print(f"\nStep 7: Verifying data...")
print(f"  Total counts per cell (first 3):")
cell_totals = np.array(adata.X.sum(axis=1)).flatten()[:3]
for i, total in enumerate(cell_totals):
    print(f"    Cell {i+1}: {total:,.0f} counts")

# Save
output_file = output_dir / "GSE132465.h5ad"
print(f"\nStep 8: Saving to {output_file}...")
adata.write(output_file)
file_size = output_file.stat().st_size / 1024 / 1024
print(f"✓ Saved! File size: {file_size:.1f} MB")

print("\n" + "="*80)
print("✅ GSE132465 conversion completed successfully!")
print("="*80)
print(f"\nSummary:")
print(f"  - Input: {len(count_matrix.columns)} cells x {len(count_matrix.index)} genes")
print(f"  - Output: {adata.shape[0]} cells x {adata.shape[1]} genes")
print(f"  - Sparsity: {sparsity:.1f}%")
print(f"  - File: {output_file}")
print(f"  - Size: {file_size:.1f} MB")
