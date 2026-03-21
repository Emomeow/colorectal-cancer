#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 1: Select top-3000 highly variable genes (batch-aware) from merged dataset.

Input : processed_data/merged_5datasets_unified.h5ad  (195,090 cells x 9,009 genes)
Output: processed_data/merged_hvg3k.h5ad              (195,090 cells x ~3,000 genes)
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import os
import sys
import io
from datetime import datetime

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

print("=" * 70)
print(f"[{datetime.now().strftime('%H:%M:%S')}] HVG Selection: 9009 -> 3000 genes")
print("=" * 70)

INPUT  = 'processed_data/merged_5datasets_unified.h5ad'
OUTPUT = 'processed_data/merged_hvg3k.h5ad'
N_TOP  = 3000

# ── Load ──────────────────────────────────────────────────────────────────────
print(f"\n[1] Loading {INPUT} ...")
adata = sc.read_h5ad(INPUT)
print(f"    Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
print(f"    Batches: {adata.obs['batch'].nunique()}  "
      f"({list(adata.obs['batch'].unique())})")

# Verify raw integer counts
x_sample = adata.X[:100, :100]
if hasattr(x_sample, 'toarray'):
    x_sample = x_sample.toarray()
is_integer = np.all(x_sample == x_sample.astype(int))
print(f"    X is integer counts: {is_integer}")
if not is_integer:
    print("    WARNING: X does not look like raw counts — HVG results may differ")

# ── Store raw counts layer ────────────────────────────────────────────────────
print(f"\n[2] Storing raw counts in adata.layers['counts'] ...")
import scipy.sparse as sp
if sp.issparse(adata.X):
    adata.layers['counts'] = adata.X.copy()
else:
    adata.layers['counts'] = sp.csr_matrix(adata.X)
print(f"    Done.")

# ── HVG selection (batch-aware, seurat_v3 on raw counts) ─────────────────────
print(f"\n[3] Selecting top-{N_TOP} HVGs (flavor='seurat_v3', batch_key='batch') ...")
print(f"    seurat_v3 uses raw counts directly (no normalization needed)")

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=N_TOP,
    batch_key='batch',
    flavor='seurat_v3',
    layer='counts',
    subset=False,          # just mark, don't subset yet
    span=0.3,              # loess smoothing span (default 0.3)
)

n_hvg = adata.var['highly_variable'].sum()
print(f"    HVGs selected: {n_hvg}")

# Summary per batch
if 'highly_variable_nbatches' in adata.var.columns:
    nbatch_counts = adata.var.loc[adata.var['highly_variable'],
                                  'highly_variable_nbatches'].value_counts().sort_index()
    print(f"\n    HVG presence across batches:")
    for nb, cnt in nbatch_counts.items():
        print(f"      present in {nb} batch(es): {cnt} genes")

# ── Subset ────────────────────────────────────────────────────────────────────
print(f"\n[4] Subsetting to {n_hvg} HVGs ...")
adata_hvg = adata[:, adata.var['highly_variable']].copy()
print(f"    Shape after subset: {adata_hvg.n_obs:,} x {adata_hvg.n_vars:,}")

# Keep only essential obs / var columns
print(f"\n    obs columns retained: {list(adata_hvg.obs.columns)}")
print(f"    var columns: {list(adata_hvg.var.columns)}")

# ── Save ──────────────────────────────────────────────────────────────────────
print(f"\n[5] Saving to {OUTPUT} ...")
adata_hvg.write_h5ad(OUTPUT, compression='gzip')
size_mb = os.path.getsize(OUTPUT) / 1024**2
print(f"    Saved: {size_mb:.1f} MB")

# ── Quick sanity check ────────────────────────────────────────────────────────
print(f"\n[Summary]")
print(f"  Input : {adata.n_obs:,} cells x {adata.n_vars:,} genes")
print(f"  Output: {adata_hvg.n_obs:,} cells x {adata_hvg.n_vars:,} genes")
print(f"  layers['counts'] present: {'counts' in adata_hvg.layers}")
print(f"\nNext step: python train_scvi_hvg3k.py")
print(f"[{datetime.now().strftime('%H:%M:%S')}] Done.")
