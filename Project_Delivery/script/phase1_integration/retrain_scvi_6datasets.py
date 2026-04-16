#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
6-Dataset UMAP Integration (正确方式)

scVI 模型仅用5个训练集训练（scvi_model_optimized），
GSE146771 作为外部验证集，通过已训练的 encoder 投影到 latent space，
不参与训练，不影响模型权重。

流程:
  1. 加载 5-dataset 整合数据（已有 X_scvi）
  2. 加载 GSE146771，TPM -> 近似计数，对齐到9009个基因空间
  3. 用 scvi_model_optimized encoder 投影 GSE146771 -> X_scvi
  4. 合并两部分 X_scvi，重新计算 UMAP（6数据集联合可视化）
  5. 保存结果和图
"""

import sys, io, os
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import scvi
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime

OUT_DIR = 'scvi_results_6datasets'
os.makedirs(OUT_DIR, exist_ok=True)

GROUPED_MAP = {
    'T cells': 'T cells', 'B cells': 'B cells',
    'Myeloids': 'Myeloids', 'NK cells': 'NK cells', 'Mast cells': 'Mast cells',
    'malignant cell': 'Malignant',
    'Epithelial cells': 'Epithelial', 'early colonocyte': 'Epithelial',
    'colonocyte': 'Epithelial', 'BEST4+ colonocyte': 'Epithelial',
    'colon goblet cell': 'Epithelial', 'secretory cell': 'Epithelial',
    'tuft cell of colon': 'Epithelial',
    'enteroendocrine cell of colon': 'Epithelial',
    'intestinal crypt stem cell of colon': 'Epithelial',
    'Stromal cells': 'Stromal',
    'Unknown': 'Unknown', 'Unspecified': 'Unknown',
}

GSE146771_TYPE_MAP = {
    'CD4 T cell': 'T cells', 'CD8 T cell': 'T cells',
    'B cell': 'B cells', 'Myeloid cell': 'Myeloids', 'ILC': 'NK cells',
}

print("=" * 70)
print(f"6-Dataset UMAP (5 trained + 1 projected)  |  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("NOTE: scVI model trained on 5 datasets only. GSE146771 = projection only.")
print("=" * 70)

# ─────────────────────────────────────────────────────────────────────────────
# 1. 加载 5-dataset 整合数据（scvi_model_optimized 的输出）
# ─────────────────────────────────────────────────────────────────────────────
print("\n[1/5] Loading 5-dataset integrated data ...")
adata5 = sc.read_h5ad('processed_data/merged_annotated.h5ad')
print(f"  Shape: {adata5.n_obs:,} cells x {adata5.n_vars:,} genes")
print(f"  X_scvi: {adata5.obsm['X_scvi'].shape}")

# 确保有 cell_type_grouped
if 'cell_type_grouped' not in adata5.obs.columns:
    adata5.obs['cell_type_grouped'] = (
        adata5.obs['unified_cell_type'].astype(str)
              .map(GROUPED_MAP).fillna('Unknown')
    )

# ─────────────────────────────────────────────────────────────────────────────
# 2. 加载 GSE146771，TPM -> 近似计数，对齐基因空间
# ─────────────────────────────────────────────────────────────────────────────
print("\n[2/5] Preparing GSE146771 for projection ...")
gse146 = sc.read_h5ad('processed_data/GSE146771_tpm.h5ad')

# TPM -> 近似 UMI 计数
raw_numi = gse146.obs['raw.nUMI'].values.astype(np.float32)
X_tpm    = gse146.X.toarray().astype(np.float32) if sp.issparse(gse146.X) else gse146.X.astype(np.float32)
X_counts = np.round(
    np.clip(np.expm1(X_tpm) * raw_numi[:, None] / 1e6, 0, None)
).astype(np.float32)

# 对齐到训练基因空间（9,009个基因）
train_genes   = adata5.var_names
common_genes  = train_genes.intersection(gse146.var_names)
missing_genes = train_genes.difference(gse146.var_names)
print(f"  Training genes : {len(train_genes):,}")
print(f"  Common genes   : {len(common_genes):,}")
print(f"  Padded with 0  : {len(missing_genes):,}")

# 构建对齐后的 AnnData（顺序必须与 train_genes 完全一致）
X_common = sp.csr_matrix(
    X_counts[:, gse146.var_names.get_indexer(common_genes)]
)
X_zero   = sp.csr_matrix((gse146.n_obs, len(missing_genes)), dtype=np.float32)

obs_shared = gse146.obs[['raw.nUMI']].copy()
obs_shared['batch']             = 'GSE132257'   # 代理 batch（最近 CRC 批次）
obs_shared['disease']           = 'CRC'
obs_shared['unified_cell_type'] = (
    gse146.obs['Global_Cluster'].astype(str).map(GSE146771_TYPE_MAP).fillna('Unknown')
)
obs_shared['cell_type_grouped'] = obs_shared['unified_cell_type']
obs_shared['real_batch']        = 'GSE146771'   # 真实来源（用于可视化）
obs_shared['unified_cell_subtype'] = gse146.obs['Sub_Cluster'].astype(str)

gse146_common = ad.AnnData(X=X_common, obs=obs_shared.copy(),
                            var=pd.DataFrame(index=common_genes))
gse146_zero   = ad.AnnData(X=X_zero,   obs=obs_shared.copy(),
                            var=pd.DataFrame(index=missing_genes))
gse146_full   = ad.concat([gse146_common, gse146_zero], axis=1)[:, train_genes].copy()
gse146_full.obs = obs_shared.copy()

print(f"  Aligned GSE146771: {gse146_full.shape}")

# ─────────────────────────────────────────────────────────────────────────────
# 3. 用 scvi_model_optimized 投影 GSE146771
# ─────────────────────────────────────────────────────────────────────────────
print("\n[3/5] Projecting GSE146771 via scVI encoder (no retraining) ...")

scvi.model.SCVI.setup_anndata(
    gse146_full,
    batch_key                  = 'batch',
    categorical_covariate_keys = ['disease'],
    layer                      = None,
)
model = scvi.model.SCVI.load('scvi_model_optimized', adata=adata5)
gse146_latent = model.get_latent_representation(adata=gse146_full)
print(f"  GSE146771 latent: {gse146_latent.shape}")

# ─────────────────────────────────────────────────────────────────────────────
# 4. 合并 X_scvi，重新计算 UMAP
# ─────────────────────────────────────────────────────────────────────────────
print("\n[4/5] Combining latents and computing joint UMAP ...")

# 构建合并用的轻量 AnnData（只需 latent embedding + obs）
obs5 = adata5.obs[['batch', 'disease', 'unified_cell_type',
                    'cell_type_grouped', 'unified_cell_subtype']].copy()
obs5['real_batch'] = obs5['batch']

obs6 = pd.DataFrame({
    'batch'               : 'GSE132257',   # 代理（模型需要）
    'disease'             : 'CRC',
    'unified_cell_type'   : obs_shared['unified_cell_type'].values,
    'cell_type_grouped'   : obs_shared['cell_type_grouped'].values,
    'unified_cell_subtype': obs_shared['unified_cell_subtype'].values,
    'real_batch'          : 'GSE146771',
}, index=gse146_full.obs_names)

obs_all   = pd.concat([obs5, obs6], axis=0)
latent_all = np.vstack([adata5.obsm['X_scvi'], gse146_latent])

combined = ad.AnnData(
    X   = sp.csr_matrix((len(obs_all), 1)),   # 占位符
    obs = obs_all,
)
combined.obsm['X_scvi'] = latent_all

sc.pp.neighbors(combined, use_rep='X_scvi', n_neighbors=15)
sc.tl.umap(combined, min_dist=0.3)
print(f"  Combined UMAP: {combined.obsm['X_umap'].shape}")
print(f"  real_batch distribution:")
for b, n in combined.obs['real_batch'].value_counts().items():
    print(f"    {b}: {n:,}")

# 保存
out_h5ad = f'{OUT_DIR}/merged_6datasets_integrated.h5ad'
combined.write_h5ad(out_h5ad, compression='gzip')
print(f"  Saved: {out_h5ad}")

# ─────────────────────────────────────────────────────────────────────────────
# 5. 生成 UMAP 图
# ─────────────────────────────────────────────────────────────────────────────
print("\n[5/5] Generating UMAP figures ...")

BATCH_COLORS = {
    'GSE132257' : '#E64B35',
    'GSE132465' : '#4DBBD5',
    'GSE125527' : '#00A087',
    'GSE150115' : '#3C5488',
    'cellxgene' : '#F39B7F',
    'GSE146771' : '#8491B4',   # 验证集，单独颜色
}
TYPE_COLORS = {
    'T cells'   : '#E64B35', 'B cells'    : '#4DBBD5',
    'Myeloids'  : '#00A087', 'NK cells'   : '#3C5488',
    'Mast cells': '#F39B7F', 'Malignant'  : '#B09C85',
    'Epithelial': '#91D1C2', 'Stromal'    : '#DC0000',
    'Unknown'   : '#AAAAAA',
}
DISEASE_COLORS = {'CRC': '#E64B35', 'UC': '#4DBBD5'}

umap = combined.obsm['X_umap']
obs  = combined.obs

def scatter(ax, col, cmap, title, s=0.3, alpha=0.5):
    for cat in sorted(obs[col].astype(str).unique()):
        mask = obs[col].astype(str) == cat
        ax.scatter(umap[mask, 0], umap[mask, 1],
                   c=cmap.get(cat, '#CCCCCC'), s=s,
                   alpha=alpha, label=cat, rasterized=True)
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xlabel('UMAP 1', fontsize=9); ax.set_ylabel('UMAP 2', fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    leg = ax.legend(markerscale=6, fontsize=7,
                    loc='upper right', framealpha=0.7,
                    bbox_to_anchor=(1.01, 1))
    handles = getattr(leg, 'legend_handles', None) or getattr(leg, 'legendHandles', [])
    for h in handles:
        h.set_alpha(1.0)

# 总览三联图
fig, axes = plt.subplots(1, 3, figsize=(26, 8))
fig.suptitle('scVI Integration: 5 Training Datasets + GSE146771 (Projected)',
             fontsize=13, fontweight='bold')
scatter(axes[0], 'real_batch',        BATCH_COLORS,   'Panel A — Dataset\n(GSE146771 = projected, not trained)')
scatter(axes[1], 'disease',           DISEASE_COLORS, 'Panel B — Disease')
scatter(axes[2], 'cell_type_grouped', TYPE_COLORS,    'Panel C — Cell Type (Grouped)')
plt.tight_layout()
out1 = f'{OUT_DIR}/umap_6datasets_overview.png'
plt.savefig(out1, dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  Saved: {out1}")

# 高清 Panel A
fig, ax = plt.subplots(figsize=(10, 8))
scatter(ax, 'real_batch', BATCH_COLORS,
        'Figure 1A — scVI Batch Integration\n(GSE146771 projected via encoder)', s=0.5)
plt.tight_layout()
out2 = f'{OUT_DIR}/figure1_panelA_batch.png'
plt.savefig(out2, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  Saved: {out2}")

# 高清细胞类型
fig, ax = plt.subplots(figsize=(10, 8))
scatter(ax, 'cell_type_grouped', TYPE_COLORS,
        'Figure 1 — Cell Type (Grouped)', s=0.5)
plt.tight_layout()
out3 = f'{OUT_DIR}/figure1_celltype_grouped.png'
plt.savefig(out3, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  Saved: {out3}")

print("\n" + "=" * 70)
print("DONE")
print(f"  Overview  : {out1}")
print(f"  Panel A   : {out2}")
print(f"  Cell Type : {out3}")
print(f"  h5ad      : {out_h5ad}")
print()
print("NOTE: scVI model = scvi_model_optimized (5 datasets only)")
print("      GSE146771 was PROJECTED, not trained")
print("=" * 70)
