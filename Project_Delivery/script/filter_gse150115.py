#!/usr/bin/env python3
"""
Filter GSE150115 with strict quality standards
"""

import anndata as ad
import scanpy as sc
import numpy as np
import os
from datetime import datetime

print(f"[{datetime.now().strftime('%H:%M:%S')}] Filtering GSE150115 with strict standards")
print("=" * 80)

# Load data
print("\nLoading GSE150115...")
adata = ad.read_h5ad('processed_data/GSE150115.h5ad')
print(f"Before filtering: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# Calculate QC metrics if not already present
if 'n_genes_by_counts' not in adata.obs.columns:
    print("\nCalculating QC metrics...")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                               log1p=False, inplace=True)

# Display pre-filtering statistics
print("\nPre-filtering QC statistics:")
print(f"  Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")
print(f"  Median counts/cell: {adata.obs['total_counts'].median():.0f}")
print(f"  Median MT%: {adata.obs['pct_counts_mt'].median():.2f}%")

# Apply STRICT filtering criteria
print("\nApplying STRICT filtering criteria:")
min_genes = 500
min_counts = 1000
max_mt = 20

print(f"  min_genes: {min_genes}")
print(f"  min_counts: {min_counts}")
print(f"  max_mt_pct: {max_mt}%")

# Create filter mask
filter_mask = (
    (adata.obs['n_genes_by_counts'] >= min_genes) &
    (adata.obs['total_counts'] >= min_counts) &
    (adata.obs['pct_counts_mt'] <= max_mt)
)

n_filtered = np.sum(~filter_mask)
pct_filtered = 100 * n_filtered / adata.n_obs

print(f"\nFiltering results:")
print(f"  Cells removed: {n_filtered:,} ({pct_filtered:.1f}%)")
print(f"  Cells retained: {np.sum(filter_mask):,} ({100-pct_filtered:.1f}%)")

# Apply filter
adata_filtered = adata[filter_mask, :].copy()

# Display post-filtering statistics
print("\nPost-filtering QC statistics:")
print(f"  Median genes/cell: {adata_filtered.obs['n_genes_by_counts'].median():.0f}")
print(f"  Median counts/cell: {adata_filtered.obs['total_counts'].median():.0f}")
print(f"  Median MT%: {adata_filtered.obs['pct_counts_mt'].median():.2f}%")

# Add dataset label
adata_filtered.obs['dataset'] = 'GSE150115'

# Save filtered data
print("\nSaving filtered data...")
output_file = 'processed_data/GSE150115_filtered.h5ad'
adata_filtered.write_h5ad(output_file, compression='gzip')

file_size_mb = os.path.getsize(output_file) / (1024**2)
print(f"Saved: {output_file}")
print(f"File size: {file_size_mb:.1f} MB")

print("\n" + "=" * 80)
print(f"[{datetime.now().strftime('%H:%M:%S')}] GSE150115 filtering complete!")
print("=" * 80)
print(f"Final: {adata_filtered.n_obs:,} cells × {adata_filtered.n_vars:,} genes")
