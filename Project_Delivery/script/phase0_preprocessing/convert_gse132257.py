"""
Convert GSE132257 to h5ad format
Test with small dataset first
"""

import pandas as pd
import scanpy as sc
import anndata as ad
from pathlib import Path
import sys

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8')

data_dir = Path("data")
output_dir = Path("processed_data")
output_dir.mkdir(exist_ok=True)

print("="*80)
print("Converting GSE132257 to h5ad format")
print("="*80)

# Load annotation
print("\nStep 1: Loading annotation...")
annot_file = data_dir / "GSE132257_processed_protocol_and_fresh_frozen_cell_annotation.txt.gz"
annotation = pd.read_csv(annot_file, compression='gzip', sep='\t', index_col=0)
print(f"✓ Loaded {len(annotation)} cells")
print(f"  Columns: {list(annotation.columns)}")
print(f"  First cell ID: {annotation.index[0]}")

# Load count matrix (genes x cells)
print("\nStep 2: Loading count matrix...")
count_file = data_dir / "GSE132257_GEO_processed_protocol_and_fresh_frozen_raw_UMI_count_matrix.txt.gz"
print(f"  File size: {count_file.stat().st_size / 1024 / 1024:.1f} MB")

# Read count matrix
count_matrix = pd.read_csv(count_file, compression='gzip', sep='\t', index_col=0)
print(f"✓ Loaded matrix: {count_matrix.shape[0]} genes x {count_matrix.shape[1]} cells")
print(f"  First cell ID in matrix: {count_matrix.columns[0]}")
print(f"  First gene: {count_matrix.index[0]}")

# Check data type
print(f"\nStep 3: Checking data type...")
print(f"  Matrix dtype: {count_matrix.dtypes[0]}")
print(f"  Sample values:\n{count_matrix.iloc[:3, :3]}")

# Convert cell IDs (dots to hyphens)
print(f"\nStep 4: Converting cell IDs (dots -> hyphens)...")
original_ids = count_matrix.columns[:3].tolist()
count_matrix.columns = count_matrix.columns.str.replace('.', '-')
converted_ids = count_matrix.columns[:3].tolist()
print(f"  Examples:")
for orig, conv in zip(original_ids, converted_ids):
    print(f"    {orig} -> {conv}")

# Match cells
print(f"\nStep 5: Matching cells...")
common_cells = list(set(count_matrix.columns) & set(annotation.index))
print(f"  Common cells: {len(common_cells)} / {len(count_matrix.columns)}")

if len(common_cells) < len(count_matrix.columns):
    print(f"  ⚠ Warning: {len(count_matrix.columns) - len(common_cells)} cells not matched")

# Create AnnData object with sparse matrix
print(f"\nStep 6: Creating AnnData object with SPARSE matrix...")
from scipy.sparse import csr_matrix

# Transpose to cells x genes
dense_matrix = count_matrix.T.loc[common_cells, :]
print(f"  Converting to sparse matrix...")
sparse_matrix = csr_matrix(dense_matrix.values)

# Calculate sparsity
n_zeros = (dense_matrix.values == 0).sum()
total = dense_matrix.values.size
sparsity = (n_zeros / total) * 100
print(f"  Sparsity: {sparsity:.1f}% zeros")

# Create AnnData
adata = ad.AnnData(X=sparse_matrix)
adata.obs_names = dense_matrix.index
adata.var_names = dense_matrix.columns
adata.obs = annotation.loc[common_cells, :]
adata.obs['dataset'] = 'GSE132257'
adata.obs['batch'] = 'GSE132257'

# Store raw counts (also sparse)
adata.layers['counts'] = adata.X.copy()

print(f"✓ Created AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes")
print(f"  X type: {type(adata.X)}")
print(f"  X dtype: {adata.X.dtype}")
print(f"  Obs columns: {list(adata.obs.columns)}")

# Verify data
print(f"\nStep 7: Verifying data...")
print(f"  Total counts per cell (first 3):")
import numpy as np
cell_totals = np.array(adata.X.sum(axis=1)).flatten()[:3]
for i, total in enumerate(cell_totals):
    print(f"    Cell {i+1}: {total:,.0f} counts")

# Save
output_file = output_dir / "GSE132257.h5ad"
print(f"\nStep 8: Saving to {output_file}...")
adata.write(output_file)
file_size = output_file.stat().st_size / 1024 / 1024
print(f"✓ Saved! File size: {file_size:.1f} MB")

print("\n" + "="*80)
print("✅ GSE132257 conversion completed successfully!")
print("="*80)
print(f"\nSummary:")
print(f"  - Input: {len(count_matrix.columns)} cells x {len(count_matrix.index)} genes")
print(f"  - Output: {adata.shape[0]} cells x {adata.shape[1]} genes")
print(f"  - File: {output_file}")
print(f"  - Size: {file_size:.1f} MB")
