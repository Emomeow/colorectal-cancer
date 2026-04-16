#!/usr/bin/env python3
"""
Convert GSE125527 to h5ad format - CORRECTED VERSION
File structure: cells × genes (rows=cells, columns=genes)
Target: 32GB memory using sparse matrices
"""

import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix, vstack
import anndata as ad
import glob
import os
import gc
from datetime import datetime

print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Starting GSE125527 conversion (CORRECTED)")
print(f"Working directory: {os.getcwd()}")

# Step 1: Load full gene list
print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Step 1: Loading gene list...")
gene_list_file = 'data/GSE125527_gene_id_rownames.csv.gz'
full_genes = pd.read_csv(gene_list_file, compression='gzip', header=None)[0].tolist()
n_genes = len(full_genes)
print(f"Total genes in reference: {n_genes:,}")

# Create gene to index mapping for O(1) lookup
gene_to_idx = {gene: idx for idx, gene in enumerate(full_genes)}
print(f"Created gene index mapping")

# Step 2: Get all sample files (ONLY UMI count files, not BCR/TCR)
print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Step 2: Finding sample files...")
sample_files = sorted(glob.glob('data/GSE125527_extracted/*UMI_table.tsv.gz'))
n_samples = len(sample_files)
print(f"Found {n_samples} UMI count sample files")

if n_samples == 0:
    raise ValueError("No sample files found!")

# Step 3: Process each sample and collect sparse matrices
print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Step 3: Processing samples...")
all_sparse_matrices = []
all_cell_ids = []
total_cells = 0

for i, sample_file in enumerate(sample_files, 1):
    sample_name = os.path.basename(sample_file).replace('_cell-gene_UMI_table.tsv.gz', '')
    print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Processing {i}/{n_samples}: {sample_name}")

    # Read sample file - cells are rows, genes are columns
    # index_col=0 sets cell_id as index
    sample_df = pd.read_csv(sample_file, compression='gzip', sep='\t', index_col=0)

    n_cells_in_sample = sample_df.shape[0]
    n_genes_in_sample = sample_df.shape[1]
    print(f"  Sample dimensions: {n_cells_in_sample} cells × {n_genes_in_sample} genes")

    # Get genes present in this sample
    sample_genes = sample_df.columns.tolist()

    # IMPORTANT: Handle gene name format differences
    # Sample files use MT.XX (dot) but reference uses MT-XX (dash)
    # Normalize gene names to match reference format
    sample_genes_normalized = []
    for gene in sample_genes:
        # Convert MT.XX to MT-XX
        if gene.startswith('MT.'):
            gene_normalized = gene.replace('MT.', 'MT-')
        else:
            gene_normalized = gene
        sample_genes_normalized.append(gene_normalized)

    # Create a full matrix aligned to reference genes (cells × full_genes)
    # Initialize with zeros in sparse format
    rows = []
    cols = []
    data = []

    # Track MT gene conversion for first sample only
    if i == 1:
        mt_genes_in_sample = [g for g in sample_genes if g.startswith('MT.')]
        if mt_genes_in_sample:
            print(f"  Found {len(mt_genes_in_sample)} MT genes with '.' format, converting to '-' format")

    for gene_idx_in_sample, gene_name_normalized in enumerate(sample_genes_normalized):
        if gene_name_normalized in gene_to_idx:
            ref_gene_idx = gene_to_idx[gene_name_normalized]  # FIXED: use normalized name
            # Get all non-zero values for this gene across all cells
            gene_values = sample_df.iloc[:, gene_idx_in_sample].values
            nonzero_mask = gene_values != 0

            if np.any(nonzero_mask):
                cell_indices = np.where(nonzero_mask)[0]
                rows.extend(cell_indices)
                cols.extend([ref_gene_idx] * len(cell_indices))
                data.extend(gene_values[nonzero_mask])

    # Create sparse matrix for this sample
    sample_sparse = csr_matrix(
        (data, (rows, cols)),
        shape=(n_cells_in_sample, n_genes),
        dtype=np.float32
    )

    sparsity = 100 * (1 - sample_sparse.nnz / (sample_sparse.shape[0] * sample_sparse.shape[1]))
    print(f"  Sparse matrix: {sample_sparse.shape[0]:,} × {sample_sparse.shape[1]:,}, "
          f"{sample_sparse.nnz:,} non-zero, {sparsity:.1f}% sparse")

    all_sparse_matrices.append(sample_sparse)

    # Add sample prefix to cell IDs to avoid conflicts
    # Extract sample ID from filename (e.g., GSM3576396_C9_R -> C9_R)
    sample_id = '_'.join(sample_name.split('_')[1:])  # Remove GSM number
    sample_cell_ids = [f"{sample_id}_{cell_id}" for cell_id in sample_df.index]
    all_cell_ids.extend(sample_cell_ids)

    total_cells += n_cells_in_sample

    # Clean up
    del sample_df, sample_sparse
    gc.collect()

    if (i % 10 == 0) or (i == n_samples):
        print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Progress: {i}/{n_samples} samples, "
              f"{total_cells:,} cells accumulated")

# Step 4: Vertically stack all samples
print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Step 4: Stacking all samples...")
print(f"Stacking {len(all_sparse_matrices)} sparse matrices...")
X_sparse = vstack(all_sparse_matrices, format='csr', dtype=np.float32)

overall_sparsity = 100 * (1 - X_sparse.nnz / (X_sparse.shape[0] * X_sparse.shape[1]))
print(f"Final matrix: {X_sparse.shape[0]:,} cells × {X_sparse.shape[1]:,} genes")
print(f"Non-zero elements: {X_sparse.nnz:,}")
print(f"Sparsity: {overall_sparsity:.2f}%")

# Clean up intermediate matrices
del all_sparse_matrices
gc.collect()

# Step 5: Create AnnData object
print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Step 5: Creating AnnData object...")
adata = ad.AnnData(
    X=X_sparse,
    obs=pd.DataFrame(index=all_cell_ids),
    var=pd.DataFrame(index=full_genes)
)

print(f"AnnData created: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# Step 6: Save to h5ad
print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Step 6: Saving to h5ad...")
output_dir = 'processed_data'
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, 'GSE125527.h5ad')

adata.write_h5ad(output_file, compression='gzip')

file_size_mb = os.path.getsize(output_file) / (1024**2)
print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] SUCCESS!")
print(f"Output: {output_file}")
print(f"File size: {file_size_mb:.1f} MB")
print(f"Final dimensions: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
print(f"Sparsity: {overall_sparsity:.2f}%")
print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Conversion complete!")
