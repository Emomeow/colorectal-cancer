#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 1 preprocessing: MT gene filtering + Scrublet doublet removal

Processes all 5 training datasets before merging:
  - GSE132257, GSE132465, GSE125527, GSE150115, cellxgene

Outputs: processed_data/{dataset}_preprocessed.h5ad
"""

import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import os
import sys
from datetime import datetime

if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

print(f"[{datetime.now().strftime('%H:%M:%S')}] Phase 1 preprocessing: MT filter + Scrublet")
print("=" * 80)
print(f"scanpy version: {sc.__version__}")

# ── Dataset configurations ─────────────────────────────────────────────────────
# Each entry describes the input file, output path, and filter thresholds.
# 'extra_filters' keys:
#   min_genes, min_counts : lower-bound cell filters (for GSE150115)
#   max_genes_mad_factor  : upper-bound = median + factor*MAD (for cellxgene)
datasets = [
    {
        'name': 'GSE132257',
        'input': 'processed_data/GSE132257.h5ad',
        'output': 'processed_data/GSE132257_preprocessed.h5ad',
        'max_mt': 10,
        # Smillie et al. annotation has 'Sample' column → per-sample Scrublet
        'scrublet_batch_candidates': ['Sample', 'sample'],
    },
    {
        'name': 'GSE132465',
        'input': 'processed_data/GSE132465.h5ad',
        'output': 'processed_data/GSE132465_preprocessed.h5ad',
        'max_mt': 10,
        # Zhang et al. annotation has 'Patient' column → per-patient Scrublet
        'scrublet_batch_candidates': ['Patient', 'Sample', 'patient', 'sample'],
    },
    {
        'name': 'GSE125527',
        'input': 'processed_data/GSE125527.h5ad',
        'output': 'processed_data/GSE125527_preprocessed.h5ad',
        'max_mt': 10,
        # Cell IDs: {sample_id}_{barcode} — sample derived from prefix
        'scrublet_batch_candidates': ['_derive_from_cell_id_gse125527'],
    },
    {
        'name': 'GSE150115',
        'input': 'processed_data/GSE150115.h5ad',
        'output': 'processed_data/GSE150115_preprocessed.h5ad',
        'max_mt': 10,
        'min_genes': 500,
        'min_counts': 1000,
        # Cell IDs: {GSM_ID}_{barcode} — GSM ID is sample
        'scrublet_batch_candidates': ['_derive_from_cell_id_gse150115'],
    },
    {
        'name': 'cellxgene',
        'input': 'processed_data/cellxgene_raw.h5ad',
        'output': 'processed_data/cellxgene_preprocessed.h5ad',
        'max_mt': 10,
        'min_genes': 200,
        'min_counts': 500,
        'max_genes_mad_factor': 3,
        # Standard CXG field; fall back to whole-dataset if absent
        'scrublet_batch_candidates': ['donor_id', 'sample_id', 'donor'],
    },
]

summary_rows = []

for cfg in datasets:
    name = cfg['name']
    print(f"\n{'=' * 60}")
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Processing: {name}")
    print(f"{'=' * 60}")

    # ── Load ──────────────────────────────────────────────────────────────────
    print(f"\nLoading {cfg['input']}...")
    adata = ad.read_h5ad(cfg['input'])
    n_cells_start = adata.n_obs
    print(f"  Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    print(f"  Obs columns: {list(adata.obs.columns[:10])}")

    # ── QC metrics ───────────────────────────────────────────────────────────
    if 'n_genes_by_counts' not in adata.obs.columns:
        print(f"\nCalculating QC metrics...")
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                                   log1p=False, inplace=True)
    else:
        # Ensure mt var column exists for display
        if 'mt' not in adata.var.columns:
            adata.var['mt'] = adata.var_names.str.startswith('MT-')

    n_mt_genes = adata.var['mt'].sum() if 'mt' in adata.var.columns else 0
    print(f"\nQC before filtering:")
    print(f"  MT genes detected: {n_mt_genes}")
    print(f"  Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")
    print(f"  Median UMI/cell:   {adata.obs['total_counts'].median():.0f}")
    print(f"  Median MT%:        {adata.obs['pct_counts_mt'].median():.2f}%")

    # ── Step 1: Cell QC filtering ─────────────────────────────────────────────
    print(f"\nStep 1: Cell QC filtering...")
    filter_mask = pd.Series(True, index=adata.obs_names)

    # MT threshold
    max_mt = cfg['max_mt']
    mt_mask = adata.obs['pct_counts_mt'] <= max_mt
    n_high_mt = (~mt_mask).sum()
    print(f"  Cells with MT% > {max_mt}%: {n_high_mt:,} ({100*n_high_mt/adata.n_obs:.1f}%) → removed")
    filter_mask &= mt_mask

    # Min genes (optional)
    if 'min_genes' in cfg:
        mg_mask = adata.obs['n_genes_by_counts'] >= cfg['min_genes']
        n_rm = (~mg_mask).sum()
        print(f"  Cells with genes < {cfg['min_genes']}: {n_rm:,} → removed")
        filter_mask &= mg_mask

    # Min counts (optional)
    if 'min_counts' in cfg:
        mc_mask = adata.obs['total_counts'] >= cfg['min_counts']
        n_rm = (~mc_mask).sum()
        print(f"  Cells with UMI < {cfg['min_counts']}: {n_rm:,} → removed")
        filter_mask &= mc_mask

    # Max genes via MAD (optional, for cellxgene)
    if 'max_genes_mad_factor' in cfg:
        factor = cfg['max_genes_mad_factor']
        med_g = adata.obs['n_genes_by_counts'].median()
        mad_g = np.median(np.abs(adata.obs['n_genes_by_counts'] - med_g))
        max_genes_thresh = med_g + factor * mad_g
        mg_max_mask = adata.obs['n_genes_by_counts'] <= max_genes_thresh
        n_rm = (~mg_max_mask).sum()
        print(f"  Cells with genes > {max_genes_thresh:.0f} (median+{factor}*MAD): {n_rm:,} → removed")
        filter_mask &= mg_max_mask

    adata = adata[filter_mask, :].copy()
    n_after_qc = adata.n_obs
    print(f"  After QC filters: {n_after_qc:,} cells (removed {n_cells_start - n_after_qc:,})")

    # ── Step 2: Scrublet doublet removal ──────────────────────────────────────
    print(f"\nStep 2: Scrublet doublet detection...")

    # Determine batch key for per-sample scrublet
    batch_key_used = None
    candidates = cfg.get('scrublet_batch_candidates', [])

    for cand in candidates:
        if cand == '_derive_from_cell_id_gse125527':
            # Cell IDs: {sample_id}_{barcode}, barcode has no underscores (10x format)
            derived = ['_'.join(c.split('_')[:-1]) for c in adata.obs_names]
            adata.obs['_scrublet_sample'] = derived
            n_s = adata.obs['_scrublet_sample'].nunique()
            examples = list(adata.obs['_scrublet_sample'].unique()[:4])
            print(f"  Derived sample IDs from cell IDs: {n_s} samples")
            print(f"  Examples: {examples}")
            batch_key_used = '_scrublet_sample'
            break

        elif cand == '_derive_from_cell_id_gse150115':
            # Cell IDs: {GSM_ID}_{barcode}
            derived = [c.split('_')[0] for c in adata.obs_names]
            adata.obs['_scrublet_sample'] = derived
            n_s = adata.obs['_scrublet_sample'].nunique()
            examples = list(adata.obs['_scrublet_sample'].unique()[:4])
            print(f"  Derived GSM IDs from cell IDs: {n_s} samples")
            print(f"  Examples: {examples}")
            batch_key_used = '_scrublet_sample'
            break

        elif cand in adata.obs.columns:
            n_batches = adata.obs[cand].nunique()
            if 1 < n_batches <= 500:
                print(f"  Using obs column '{cand}' as batch key: {n_batches} batches")
                batch_key_used = cand
                break

    if batch_key_used is None:
        print(f"  No batch column found — running Scrublet on entire dataset")

    try:
        if batch_key_used:
            sc.pp.scrublet(adata, batch_key=batch_key_used, verbose=False)
        else:
            sc.pp.scrublet(adata, verbose=False)

        n_doublets = adata.obs['predicted_doublet'].sum()
        pct_doublets = 100 * n_doublets / adata.n_obs
        print(f"  Predicted doublets: {n_doublets:,} ({pct_doublets:.1f}%)")
        print(f"  Median doublet score: {adata.obs['doublet_score'].median():.4f}")

        adata = adata[~adata.obs['predicted_doublet'], :].copy()
        # Remove the temporary sample column if added
        if '_scrublet_sample' in adata.obs.columns:
            adata.obs.drop(columns=['_scrublet_sample'], inplace=True)

    except Exception as e:
        print(f"  WARNING: Scrublet failed: {e}")
        print(f"  Skipping doublet removal for {name}")
        n_doublets = 0

    n_after_scrublet = adata.n_obs
    print(f"  After doublet removal: {n_after_scrublet:,} cells")

    # ── Summary ───────────────────────────────────────────────────────────────
    pct_retained = 100 * n_after_scrublet / n_cells_start
    print(f"\nSummary for {name}:")
    print(f"  Start:           {n_cells_start:,} cells")
    print(f"  After QC filter: {n_after_qc:,} cells  (removed {n_cells_start - n_after_qc:,})")
    print(f"  After Scrublet:  {n_after_scrublet:,} cells  (removed {n_after_qc - n_after_scrublet:,})")
    print(f"  Final retained:  {pct_retained:.1f}%")

    # ── Save ──────────────────────────────────────────────────────────────────
    os.makedirs('processed_data', exist_ok=True)
    print(f"\nSaving to {cfg['output']}...")
    adata.write_h5ad(cfg['output'], compression='gzip')
    file_size_mb = os.path.getsize(cfg['output']) / (1024 ** 2)
    print(f"  Saved: {file_size_mb:.1f} MB")

    summary_rows.append({
        'dataset': name,
        'cells_start': n_cells_start,
        'cells_after_qc': n_after_qc,
        'n_doublets_removed': n_after_qc - n_after_scrublet,
        'cells_final': n_after_scrublet,
        'pct_retained': round(pct_retained, 1),
        'output': cfg['output'],
    })

# ── Final summary ─────────────────────────────────────────────────────────────
print(f"\n{'=' * 80}")
print(f"[{datetime.now().strftime('%H:%M:%S')}] Phase 1 preprocessing complete!")
print(f"{'=' * 80}")

summary_df = pd.DataFrame(summary_rows)
print(f"\nFinal cell counts per dataset:")
for _, row in summary_df.iterrows():
    print(f"  {row['dataset']:20s}: {row['cells_start']:>7,} → {row['cells_final']:>7,} cells "
          f"({row['pct_retained']:.1f}% retained)")

total_final = summary_df['cells_final'].sum()
print(f"\n  {'Total':20s}: {'':>8} → {total_final:>7,} cells")

summary_df.to_csv('processed_data/preprocessing_summary.csv', index=False)
print(f"\nSaved summary: processed_data/preprocessing_summary.csv")
print(f"\nNext step: run remerge_with_unified_annotations.py")
