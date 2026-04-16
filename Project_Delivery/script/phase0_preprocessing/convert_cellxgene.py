#!/usr/bin/env python3
"""
Convert cellxgene validation dataset to standardized h5ad format
Extract raw UMI counts from .raw attribute
"""

import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
from datetime import datetime
import os

print(f"[{datetime.now().strftime('%H:%M:%S')}] Processing cellxgene validation dataset")
print("=" * 80)

# Load cellxgene data
print("\nLoading cellxgene.h5ad...")
adata = ad.read_h5ad('data/cellxgene.h5ad')

print(f"Initial dimensions: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
print(f"Data type: {type(adata.X)}")
print(f"Data range: [{adata.X.min():.2f}, {adata.X.max():.2f}]")

# Check if data is normalized
if adata.X.max() < 100:
    print("\nWARNING: Data appears to be log-normalized (max < 100)")
    print("   Extracting raw counts from .raw attribute...")

    if adata.raw is not None:
        # Extract raw counts
        raw_X = adata.raw.X
        print(f"\nRaw counts found!")
        print(f"  Type: {type(raw_X)}")
        print(f"  Range: [{raw_X.min():.0f}, {raw_X.max():.0f}]")

        # Check if integer (UMI counts)
        if hasattr(raw_X, 'data'):
            is_integer = np.allclose(raw_X.data, np.round(raw_X.data))
        else:
            is_integer = np.allclose(raw_X, np.round(raw_X))

        if is_integer:
            print(f"  OK: Data contains integer UMI counts")
        else:
            print(f"  WARNING:  WARNING: Data is not integer counts!")

        # Create new AnnData with raw counts
        adata_raw = ad.AnnData(
            X=raw_X,
            obs=adata.obs.copy(),
            var=adata.raw.var.copy()
        )

        print(f"\nCreated new AnnData with raw counts:")
        print(f"  Dimensions: {adata_raw.n_obs:,} cells × {adata_raw.n_vars:,} genes")
        print(f"  Gene name format: {adata_raw.var_names[0]} (ENSEMBL IDs)")

        # Convert ENSEMBL IDs to gene symbols
        print("\nConverting ENSEMBL IDs to gene symbols...")
        if 'feature_name' in adata_raw.var.columns:
            print("  Found 'feature_name' column with gene symbols")
            gene_symbols = adata_raw.var['feature_name'].values

            # Remove duplicates by keeping first occurrence
            unique_mask = ~pd.Series(gene_symbols).duplicated(keep='first').values
            n_duplicates = np.sum(~unique_mask)

            if n_duplicates > 0:
                print(f"  Removing {n_duplicates} duplicate gene symbols (keeping first)")
                adata_raw = adata_raw[:, unique_mask].copy()

            # Set gene symbols as var_names
            adata_raw.var_names = adata_raw.var['feature_name'].values
            adata_raw.var_names_make_unique()

            print(f"  After conversion: {adata_raw.n_vars:,} genes")
            print(f"  New gene name format: {adata_raw.var_names[0]}")
        else:
            print("  WARNING: No gene symbol column found, keeping ENSEMBL IDs")
            print("  This may affect MT gene detection!")

        # Replace adata with raw version
        adata = adata_raw
    else:
        print("\nERROR: ERROR: No .raw attribute found!")
        raise ValueError("Cannot extract raw counts from cellxgene data")

# Add dataset label
adata.obs['dataset'] = 'cellxgene'

# Calculate QC metrics
print("\nCalculating QC metrics...")
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                           log1p=False, inplace=True)

# Display QC statistics
print("\nQC Statistics:")
print(f"  Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")
print(f"  Median UMI/cell: {adata.obs['total_counts'].median():.0f}")
print(f"  Median MT%: {adata.obs['pct_counts_mt'].median():.2f}%")

# Check cell type distribution
if 'cell_type' in adata.obs.columns:
    print(f"\nCell type distribution (top 10):")
    ct_counts = adata.obs['cell_type'].value_counts()
    for ct, count in ct_counts.head(10).items():
        pct = 100 * count / adata.n_obs
        print(f"  {ct}: {count:,} ({pct:.1f}%)")

# Check sample information
if 'Sample Type' in adata.obs.columns:
    print(f"\nSample type distribution:")
    st_counts = adata.obs['Sample Type'].value_counts()
    for st, count in st_counts.items():
        pct = 100 * count / adata.n_obs
        print(f"  {st}: {count:,} ({pct:.1f}%)")

# Check sparsity
sparsity = 100 * (1 - adata.X.nnz / (adata.n_obs * adata.n_vars))
print(f"\nSparsity: {sparsity:.2f}%")

# Save to processed_data
print("\nSaving standardized h5ad file...")
os.makedirs('processed_data', exist_ok=True)
output_file = 'processed_data/cellxgene_raw.h5ad'

adata.write_h5ad(output_file, compression='gzip')

file_size_mb = os.path.getsize(output_file) / (1024**2)
print(f"\nOK: Saved: {output_file}")
print(f"  File size: {file_size_mb:.1f} MB")
print(f"  Dimensions: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

print("\n" + "=" * 80)
print(f"[{datetime.now().strftime('%H:%M:%S')}] cellxgene conversion complete!")
print("=" * 80)

print("\nWARNING:  IMPORTANT NOTE:")
print("This is a VALIDATION dataset and should NOT be used for scVI training because:")
print("1. Different tissue source (colorectal cancer)")
print("2. Different cell type composition (mostly epithelial/tumor cells)")
print("3. Needs to remain independent for model evaluation")
