#!/usr/bin/env python3
"""
Convert GSE146771 validation dataset to h5ad format
NOTE: This dataset contains TPM-normalized data, NOT raw UMI counts
"""

import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import anndata as ad
import scanpy as sc
from datetime import datetime
import os

print(f"[{datetime.now().strftime('%H:%M:%S')}] Processing GSE146771 validation dataset")
print("=" * 80)

print("\nWARNING:  CRITICAL NOTE:")
print("This dataset contains TPM (Transcripts Per Million) normalized data,")
print("NOT raw UMI counts. This makes it unsuitable for scVI training,")
print("but useful for validation/testing of trained models.")
print("=" * 80)

# Step 1: Load metadata
print(f"\n[{datetime.now().strftime('%H:%M:%S')}] Step 1: Loading metadata...")
metadata = pd.read_csv('data/GSE146771_CRC.Leukocyte.10x.Metadata.txt.gz',
                       sep='\t', compression='gzip')

print(f"Metadata loaded: {len(metadata):,} cells")
print(f"Columns: {list(metadata.columns)}")

# Step 2: Load TPM matrix (genes × cells)
print(f"\n[{datetime.now().strftime('%H:%M:%S')}] Step 2: Loading TPM matrix...")
print("This may take a few minutes (434 MB compressed file)...")

# Read TPM matrix - first column is gene names, rest are cells
# NOTE: File is space-delimited, not tab-delimited
tpm_df = pd.read_csv('data/GSE146771_CRC.Leukocyte.10x.TPM.txt.gz',
                     sep='\\s+', compression='gzip', index_col=0)

print(f"TPM matrix loaded: {tpm_df.shape[0]:,} genes × {tpm_df.shape[1]:,} cells")

# Verify dimensions match metadata
assert tpm_df.shape[1] == len(metadata), \
    f"Mismatch: TPM has {tpm_df.shape[1]} cells but metadata has {len(metadata)}"

# Step 3: Transpose to cells × genes
print(f"\n[{datetime.now().strftime('%H:%M:%S')}] Step 3: Transposing matrix...")
tpm_df_t = tpm_df.T

# Convert to sparse matrix (TPM values can be sparse too)
print(f"Converting to sparse matrix...")
X_sparse = csr_matrix(tpm_df_t.values, dtype=np.float32)

sparsity = 100 * (1 - X_sparse.nnz / (X_sparse.shape[0] * X_sparse.shape[1]))
print(f"Sparse matrix: {X_sparse.shape[0]:,} cells × {X_sparse.shape[1]:,} genes")
print(f"Sparsity: {sparsity:.2f}%")
print(f"Non-zero elements: {X_sparse.nnz:,}")

# Step 4: Create AnnData object
print(f"\n[{datetime.now().strftime('%H:%M:%S')}] Step 4: Creating AnnData...")

# Prepare obs (cell metadata)
obs_df = metadata.copy()
obs_df.index = tpm_df_t.index
obs_df['dataset'] = 'GSE146771'

# Prepare var (gene metadata)
var_df = pd.DataFrame(index=tpm_df.index)

# Create AnnData
adata = ad.AnnData(
    X=X_sparse,
    obs=obs_df,
    var=var_df
)

print(f"AnnData created: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# Step 5: Add QC metrics from metadata
print(f"\n[{datetime.now().strftime('%H:%M:%S')}] Step 5: Processing QC metrics...")

# The metadata already has raw UMI and gene counts
if 'raw.nUMI' in adata.obs.columns:
    print("\nOriginal raw QC metrics (before normalization):")
    print(f"  Median UMI (raw): {adata.obs['raw.nUMI'].median():.0f}")
    print(f"  Median genes (raw): {adata.obs['raw.nGene'].median():.0f}")

if 'filter.nUMI' in adata.obs.columns:
    print(f"  Median UMI (filtered): {adata.obs['filter.nUMI'].median():.0f}")
    print(f"  Median genes (filtered): {adata.obs['filter.nGene'].median():.0f}")

# Calculate metrics on TPM data (for reference)
print("\nCalculating metrics on TPM-normalized data...")
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                           log1p=False, inplace=True)

print(f"  Median genes/cell (TPM): {adata.obs['n_genes_by_counts'].median():.0f}")
print(f"  Median TPM sum/cell: {adata.obs['total_counts'].median():.2f}")
print(f"  Median MT% (TPM): {adata.obs['pct_counts_mt'].median():.2f}%")

# Step 6: Cell type distribution
print(f"\n[{datetime.now().strftime('%H:%M:%S')}] Step 6: Cell type analysis...")

if 'Global_Cluster' in adata.obs.columns:
    print("\nGlobal cell type distribution:")
    gc_counts = adata.obs['Global_Cluster'].value_counts()
    for ct, count in gc_counts.items():
        pct = 100 * count / adata.n_obs
        print(f"  {ct}: {count:,} ({pct:.1f}%)")

if 'Tissue' in adata.obs.columns:
    print("\nTissue type distribution:")
    tissue_counts = adata.obs['Tissue'].value_counts()
    for tissue, count in tissue_counts.items():
        pct = 100 * count / adata.n_obs
        print(f"  {tissue}: {count:,} ({pct:.1f}%)")

# Step 7: Save
print(f"\n[{datetime.now().strftime('%H:%M:%S')}] Step 7: Saving to h5ad...")
os.makedirs('processed_data', exist_ok=True)
output_file = 'processed_data/GSE146771_tpm.h5ad'

adata.write_h5ad(output_file, compression='gzip')

file_size_mb = os.path.getsize(output_file) / (1024**2)
print(f"\nOK: Saved: {output_file}")
print(f"  File size: {file_size_mb:.1f} MB")
print(f"  Dimensions: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

print("\n" + "=" * 80)
print(f"[{datetime.now().strftime('%H:%M:%S')}] GSE146771 conversion complete!")
print("=" * 80)

print("\nWARNING:  IMPORTANT LIMITATIONS:")
print("1. Data is TPM-NORMALIZED, not raw UMI counts")
print("2. Cannot be used for scVI training (requires raw counts)")
print("3. Suitable for validation/testing only")
print("4. Original raw UMI counts are recorded in metadata but not in matrix")
