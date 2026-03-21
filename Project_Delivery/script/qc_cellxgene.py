#!/usr/bin/env python3
"""
Quality Control analysis for cellxgene dataset
"""

import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import os

print(f"[{datetime.now().strftime('%H:%M:%S')}] QC Analysis for cellxgene")
print("=" * 80)

# Load data
print("\nLoading cellxgene_raw.h5ad...")
adata = ad.read_h5ad('processed_data/cellxgene_raw.h5ad')
print(f"Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# Check if QC metrics exist
if 'n_genes_by_counts' not in adata.obs.columns:
    print("\nCalculating QC metrics...")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                               log1p=False, inplace=True)
else:
    print("\nQC metrics already present")

# Display statistics
print("\n" + "=" * 80)
print("QC Statistics")
print("=" * 80)

print(f"\nGenes per cell:")
print(f"  Median: {adata.obs['n_genes_by_counts'].median():.0f}")
print(f"  Mean: {adata.obs['n_genes_by_counts'].mean():.0f}")
print(f"  Min: {adata.obs['n_genes_by_counts'].min():.0f}")
print(f"  Max: {adata.obs['n_genes_by_counts'].max():.0f}")

print(f"\nUMI counts per cell:")
print(f"  Median: {adata.obs['total_counts'].median():.0f}")
print(f"  Mean: {adata.obs['total_counts'].mean():.0f}")
print(f"  Min: {adata.obs['total_counts'].min():.0f}")
print(f"  Max: {adata.obs['total_counts'].max():.0f}")

print(f"\nMitochondrial percentage:")
print(f"  Median: {adata.obs['pct_counts_mt'].median():.2f}%")
print(f"  Mean: {adata.obs['pct_counts_mt'].mean():.2f}%")
print(f"  Min: {adata.obs['pct_counts_mt'].min():.2f}%")
print(f"  Max: {adata.obs['pct_counts_mt'].max():.2f}%")

# Check MT genes
n_mt = adata.var_names.str.startswith('MT-').sum()
print(f"\nMitochondrial genes: {n_mt}")
if n_mt > 0:
    mt_genes = adata.var_names[adata.var_names.str.startswith('MT-')].tolist()
    print(f"  MT genes: {mt_genes}")

# Calculate suggested thresholds using MAD
print("\n" + "=" * 80)
print("Suggested Filtering Thresholds (MAD-based)")
print("=" * 80)

# Standard thresholds
min_genes = 200
min_counts = 500
max_mt = 20

# Calculate max_genes using MAD
median_genes = adata.obs['n_genes_by_counts'].median()
mad_genes = np.median(np.abs(adata.obs['n_genes_by_counts'] - median_genes))
max_genes = median_genes + 3 * mad_genes

print(f"\nSuggested thresholds:")
print(f"  min_genes: {min_genes}")
print(f"  max_genes: {max_genes:.0f} (median + 3*MAD)")
print(f"  min_counts: {min_counts}")
print(f"  max_mt_pct: {max_mt}%")

# Calculate how many cells would be filtered
filter_mask = (
    (adata.obs['n_genes_by_counts'] >= min_genes) &
    (adata.obs['n_genes_by_counts'] <= max_genes) &
    (adata.obs['total_counts'] >= min_counts) &
    (adata.obs['pct_counts_mt'] <= max_mt)
)

n_filtered = np.sum(~filter_mask)
pct_filtered = 100 * n_filtered / adata.n_obs

print(f"\nWith these thresholds:")
print(f"  Cells to filter: {n_filtered:,} ({pct_filtered:.1f}%)")
print(f"  Cells to keep: {np.sum(filter_mask):,} ({100-pct_filtered:.1f}%)")

# Create QC plots
print("\n" + "=" * 80)
print("Generating QC Plots")
print("=" * 80)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('cellxgene QC Metrics', fontsize=16, y=0.995)

# 1. Genes per cell
ax = axes[0, 0]
ax.hist(adata.obs['n_genes_by_counts'], bins=100, color='steelblue', edgecolor='black')
ax.axvline(min_genes, color='red', linestyle='--', linewidth=2, label=f'min_genes={min_genes}')
ax.axvline(max_genes, color='orange', linestyle='--', linewidth=2, label=f'max_genes={max_genes:.0f}')
ax.axvline(median_genes, color='green', linestyle='--', linewidth=2, label=f'median={median_genes:.0f}')
ax.set_xlabel('Number of genes', fontsize=12)
ax.set_ylabel('Number of cells', fontsize=12)
ax.set_title('Genes per cell', fontsize=14)
ax.legend()
ax.grid(alpha=0.3)

# 2. UMI counts per cell
ax = axes[0, 1]
ax.hist(adata.obs['total_counts'], bins=100, color='steelblue', edgecolor='black')
ax.axvline(min_counts, color='red', linestyle='--', linewidth=2, label=f'min_counts={min_counts}')
ax.axvline(adata.obs['total_counts'].median(), color='green', linestyle='--', linewidth=2,
           label=f"median={adata.obs['total_counts'].median():.0f}")
ax.set_xlabel('Total UMI counts', fontsize=12)
ax.set_ylabel('Number of cells', fontsize=12)
ax.set_title('UMI counts per cell', fontsize=14)
ax.legend()
ax.grid(alpha=0.3)

# 3. Mitochondrial percentage
ax = axes[1, 0]
ax.hist(adata.obs['pct_counts_mt'], bins=100, color='steelblue', edgecolor='black')
ax.axvline(max_mt, color='red', linestyle='--', linewidth=2, label=f'max_mt={max_mt}%')
ax.axvline(adata.obs['pct_counts_mt'].median(), color='green', linestyle='--', linewidth=2,
           label=f"median={adata.obs['pct_counts_mt'].median():.2f}%")
ax.set_xlabel('Mitochondrial %', fontsize=12)
ax.set_ylabel('Number of cells', fontsize=12)
ax.set_title('Mitochondrial percentage', fontsize=14)
ax.legend()
ax.grid(alpha=0.3)

# 4. Scatter: genes vs UMI, colored by MT%
ax = axes[1, 1]
scatter = ax.scatter(adata.obs['total_counts'],
                    adata.obs['n_genes_by_counts'],
                    c=adata.obs['pct_counts_mt'],
                    cmap='viridis',
                    s=1,
                    alpha=0.5)
ax.set_xlabel('Total UMI counts', fontsize=12)
ax.set_ylabel('Number of genes', fontsize=12)
ax.set_title('Genes vs UMI (colored by MT%)', fontsize=14)
ax.set_xscale('log')
ax.set_yscale('log')
cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('MT %', fontsize=10)
ax.grid(alpha=0.3)

plt.tight_layout()

# Save plot
os.makedirs('qc_results', exist_ok=True)
output_plot = 'qc_results/cellxgene_qc_metrics.png'
plt.savefig(output_plot, dpi=150, bbox_inches='tight')
print(f"\nSaved QC plot: {output_plot}")
plt.close()

# Create QC report CSV
print("\n" + "=" * 80)
print("Generating QC Report")
print("=" * 80)

qc_summary = {
    'dataset': 'cellxgene',
    'n_cells_before': adata.n_obs,
    'n_genes': adata.n_vars,
    'n_mt_genes': n_mt,
    'median_genes_per_cell': adata.obs['n_genes_by_counts'].median(),
    'median_counts_per_cell': adata.obs['total_counts'].median(),
    'median_mt_pct': adata.obs['pct_counts_mt'].median(),
    'suggested_min_genes': min_genes,
    'suggested_max_genes': max_genes,
    'suggested_min_counts': min_counts,
    'suggested_max_mt': max_mt,
    'n_cells_to_filter': n_filtered,
    'pct_cells_to_filter': pct_filtered
}

qc_df = pd.DataFrame([qc_summary])
output_csv = 'qc_results/cellxgene_qc_report.csv'
qc_df.to_csv(output_csv, index=False)
print(f"\nSaved QC report: {output_csv}")

# Display summary
print("\n" + "=" * 80)
print("QC Report Summary")
print("=" * 80)
print(qc_df.to_string(index=False))

print("\n" + "=" * 80)
print(f"[{datetime.now().strftime('%H:%M:%S')}] QC Analysis Complete!")
print("=" * 80)
print(f"\nOutputs:")
print(f"  1. {output_plot}")
print(f"  2. {output_csv}")
