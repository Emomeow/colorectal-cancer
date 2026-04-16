#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Remerge datasets with unified cell type annotations
"""

import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
from datetime import datetime
import os
import sys

# Set stdout encoding to UTF-8
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

print(f"[{datetime.now().strftime('%H:%M:%S')}] Remerging datasets with unified annotations")
print("=" * 80)

# Load the 5 datasets (preprocessed: MT filtered + Scrublet doublet removed)
print("\n[Step 1] Loading 5 preprocessed datasets...")
training_datasets = {
    'GSE132257': 'processed_data/GSE132257_preprocessed.h5ad',
    'GSE132465': 'processed_data/GSE132465_preprocessed.h5ad',
    'GSE125527': 'processed_data/GSE125527_preprocessed.h5ad',
    'GSE150115_filtered': 'processed_data/GSE150115_preprocessed.h5ad',
    'cellxgene': 'processed_data/cellxgene_preprocessed.h5ad'
}

adatas = {}
for name, path in training_datasets.items():
    print(f"\nLoading {name}...")
    adata = ad.read_h5ad(path)

    # 基本信息
    print(f"  维度: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # 检查并统一细胞类型列
    print(f"  可用列: {list(adata.obs.columns)[:10]}...")

    # 统一细胞类型注释到unified_cell_type列
    if 'Cell_type' in adata.obs.columns:
        adata.obs['unified_cell_type'] = adata.obs['Cell_type']
        print(f"  使用Cell_type列 -> unified_cell_type")
    elif 'cell_type' in adata.obs.columns:
        adata.obs['unified_cell_type'] = adata.obs['cell_type']
        print(f"  使用cell_type列 -> unified_cell_type")
    else:
        adata.obs['unified_cell_type'] = 'Unknown'
        print(f"  无细胞类型注释 -> 标记为Unknown")

    # 统一细胞亚型到unified_cell_subtype列
    if 'Cell_subtype' in adata.obs.columns:
        adata.obs['unified_cell_subtype'] = adata.obs['Cell_subtype']
    else:
        adata.obs['unified_cell_subtype'] = 'Unknown'

    # 添加疾病信息
    if name in ['GSE125527', 'GSE150115_filtered']:
        adata.obs['disease'] = 'UC'
    else:
        adata.obs['disease'] = 'CRC'

    # 确保dataset标签
    if 'dataset' not in adata.obs.columns:
        adata.obs['dataset'] = name.replace('_filtered', '')

    # 计算QC metrics if needed
    if 'n_genes_by_counts' not in adata.obs.columns:
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                                   log1p=False, inplace=True)

    # 显示细胞类型分布
    if adata.obs['unified_cell_type'].notna().sum() > 0:
        ct_dist = adata.obs['unified_cell_type'].value_counts()
        if len(ct_dist) > 0 and ct_dist.iloc[0] != 'Unknown':
            print(f"  细胞类型分布 (top 5):")
            for ct, count in ct_dist.head(5).items():
                if ct != 'Unknown':
                    pct = 100 * count / adata.n_obs
                    print(f"    {ct}: {count:,} ({pct:.1f}%)")

    adatas[name] = adata

# cellxgene QC filtering is now done in preprocess_phase1.py (MT + Scrublet)
print(f"\n[Step 2] cellxgene已在preprocess_phase1.py中完成QC过滤，跳过")

# Find common genes
print(f"\n[Step 3] 寻找共同基因...")
gene_sets = {name: set(adata.var_names) for name, adata in adatas.items()}

print(f"\n各数据集基因数:")
for name, genes in gene_sets.items():
    print(f"  {name}: {len(genes):,} genes")

common_genes = set.intersection(*gene_sets.values())
common_genes_sorted = sorted(common_genes)

print(f"\n共同基因: {len(common_genes_sorted):,}")

# Check MT genes
mt_common = [g for g in common_genes_sorted if g.startswith('MT-')]
print(f"  包含MT基因: {len(mt_common)} ({', '.join(mt_common)})")

# Subset to common genes
print(f"\n[Step 4] 子集到共同基因...")
subset_adatas = []
dataset_order = ['GSE132257', 'GSE132465', 'GSE125527', 'GSE150115_filtered', 'cellxgene']

for name in dataset_order:
    adata = adatas[name]
    print(f"\n{name}:")
    print(f"  子集前: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    adata_subset = adata[:, common_genes_sorted].copy()
    print(f"  子集后: {adata_subset.n_obs:,} cells × {adata_subset.n_vars:,} genes")

    subset_adatas.append(adata_subset)

# Concatenate
print(f"\n[Step 5] 合并所有数据集...")
batch_keys = ['GSE132257', 'GSE132465', 'GSE125527', 'GSE150115', 'cellxgene']

adata_merged = ad.concat(subset_adatas, join='outer', label='batch',
                         keys=batch_keys, index_unique='_')

print(f"\n合并后数据集:")
print(f"  总细胞数: {adata_merged.n_obs:,}")
print(f"  总基因数: {adata_merged.n_vars:,}")
print(f"  批次数: {len(adata_merged.obs['batch'].unique())}")

# 统计细胞类型分布
print(f"\n[Step 6] 统计细胞类型分布...")

print(f"\n按批次统计:")
batch_counts = adata_merged.obs['batch'].value_counts().sort_index()
for batch, count in batch_counts.items():
    pct = 100 * count / adata_merged.n_obs
    print(f"  {batch}: {count:,} cells ({pct:.1f}%)")

print(f"\n按疾病统计:")
disease_counts = adata_merged.obs['disease'].value_counts()
for disease, count in disease_counts.items():
    pct = 100 * count / adata_merged.n_obs
    print(f"  {disease}: {count:,} cells ({pct:.1f}%)")

print(f"\n按细胞类型统计 (top 20):")
ct_counts = adata_merged.obs['unified_cell_type'].value_counts()
for ct, count in ct_counts.head(20).items():
    pct = 100 * count / adata_merged.n_obs
    print(f"  {ct}: {count:,} ({pct:.1f}%)")

# Gene filtering
print(f"\n[Step 7] 基因过滤...")
print(f"过滤前: {adata_merged.n_vars:,} genes")
sc.pp.filter_genes(adata_merged, min_cells=10)
print(f"过滤后 (min_cells=10): {adata_merged.n_vars:,} genes")

# Calculate overall QC statistics
print(f"\n[Step 8] 计算整体QC统计...")
sparsity = 100 * (1 - adata_merged.X.nnz / (adata_merged.n_obs * adata_merged.n_vars))
print(f"稀疏度: {sparsity:.2f}%")

print(f"\nQC统计 (整体):")
print(f"  中位基因数/细胞: {adata_merged.obs['n_genes_by_counts'].median():.0f}")
print(f"  中位UMI数/细胞: {adata_merged.obs['total_counts'].median():.0f}")
print(f"  中位MT%: {adata_merged.obs['pct_counts_mt'].median():.2f}%")

# Save
print(f"\n[Step 9] 保存统一后的合并数据集...")
os.makedirs('processed_data', exist_ok=True)
output_file = 'processed_data/merged_5datasets_unified.h5ad'

adata_merged.write_h5ad(output_file, compression='gzip')

file_size_mb = os.path.getsize(output_file) / (1024**2)
print(f"\n保存成功: {output_file}")
print(f"文件大小: {file_size_mb:.1f} MB ({file_size_mb/1024:.2f} GB)")

# Generate detailed summary
print(f"\n[Step 10] 生成详细摘要...")

# Per-batch statistics
batch_stats = []
for batch in batch_keys:
    batch_data = adata_merged[adata_merged.obs['batch'] == batch]

    # 细胞类型统计
    ct_dist = batch_data.obs['unified_cell_type'].value_counts()
    top_cell_type = ct_dist.index[0] if len(ct_dist) > 0 else 'Unknown'
    n_cell_types = (ct_dist != 'Unknown').sum()

    batch_stats.append({
        'batch': batch,
        'n_cells': batch_data.n_obs,
        'pct_cells': 100 * batch_data.n_obs / adata_merged.n_obs,
        'disease': batch_data.obs['disease'].iloc[0],
        'n_cell_types': n_cell_types,
        'top_cell_type': top_cell_type,
        'median_genes': batch_data.obs['n_genes_by_counts'].median(),
        'median_umi': batch_data.obs['total_counts'].median(),
        'median_mt_pct': batch_data.obs['pct_counts_mt'].median()
    })

batch_stats_df = pd.DataFrame(batch_stats)
batch_stats_df.to_csv('processed_data/merged_unified_batch_stats.csv', index=False)
print(f"保存批次统计: processed_data/merged_unified_batch_stats.csv")

# Cell type distribution across batches
ct_batch_df = pd.crosstab(
    adata_merged.obs['unified_cell_type'],
    adata_merged.obs['batch']
)
ct_batch_df.to_csv('processed_data/merged_unified_celltype_by_batch.csv')
print(f"保存细胞类型×批次交叉表: processed_data/merged_unified_celltype_by_batch.csv")

# Overall summary
summary = {
    'total_cells': adata_merged.n_obs,
    'total_genes': adata_merged.n_vars,
    'n_batches': len(batch_keys),
    'n_diseases': len(adata_merged.obs['disease'].unique()),
    'n_cell_types': (adata_merged.obs['unified_cell_type'] != 'Unknown').sum(),
    'n_cell_types_unique': len([ct for ct in adata_merged.obs['unified_cell_type'].unique() if ct != 'Unknown']),
    'sparsity_pct': sparsity,
    'file_size_mb': file_size_mb,
    'median_genes_per_cell': adata_merged.obs['n_genes_by_counts'].median(),
    'median_umi_per_cell': adata_merged.obs['total_counts'].median(),
    'median_mt_pct': adata_merged.obs['pct_counts_mt'].median(),
    'n_mt_genes': len(mt_common),
    'CRC_cells': (adata_merged.obs['disease'] == 'CRC').sum(),
    'UC_cells': (adata_merged.obs['disease'] == 'UC').sum()
}

summary_df = pd.DataFrame([summary])
summary_df.to_csv('processed_data/merged_unified_summary.csv', index=False)
print(f"保存总体摘要: processed_data/merged_unified_summary.csv")

print("\n" + "=" * 80)
print(f"[{datetime.now().strftime('%H:%M:%S')}] 数据合并完成！")
print("=" * 80)
print(f"\n最终统一数据集: {output_file}")
print(f"维度: {adata_merged.n_obs:,} cells × {adata_merged.n_vars:,} genes")
print(f"统一列:")
print(f"  - unified_cell_type: 统一的细胞类型注释")
print(f"  - unified_cell_subtype: 统一的细胞亚型注释")
print(f"  - disease: CRC或UC")
print(f"  - batch: 5个批次")
print(f"\n准备好用于scVI训练!")
