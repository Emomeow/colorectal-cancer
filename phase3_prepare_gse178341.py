#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 3 — GSE178341 数据预处理
发现集：MMRp CRC 肿瘤细胞（T 样本）
输出：processed_data/GSE178341_MMRp.h5ad

步骤：
1. 从 H5 文件加载全量表达矩阵
2. 合并 metatables（患者信息/MMR状态）和 cluster 注释
3. 筛选 MMRp + 肿瘤（SPECIMEN_TYPE == 'T'）
4. QC：过滤低质量细胞
5. 标准化 + log1p
6. 高变基因选择 + PCA + UMAP（快速预览）
7. 打 Pre-driver 基因模块分（TAM / CD8_Tex signature from Phase 2）
8. 保存
"""

import scanpy as sc
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, csc_matrix
import h5py, gzip, io, sys, os, warnings
warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

sc.settings.verbosity = 1

OUT_DIR = 'processed_data'
os.makedirs(OUT_DIR, exist_ok=True)
OUT_H5AD = f'{OUT_DIR}/GSE178341_MMRp.h5ad'

# ── Phase 2 Pre-driver signature genes（与 GSE144735 打分一致）──────
TAM_MODULE = ['FTH1', 'SOD2', 'SOCS3', 'OSM', 'HMOX1', 'FCGR1A',
              'LGALS3', 'CD68', 'HLA-DPB1', 'HLA-DPA1', 'RGS1',
              'CTSB', 'HLA-DRA', 'HLA-DRB1', 'GPR183', 'APOE', 'SPP1']

TEX_MODULE  = ['HAVCR2', 'ENTPD1', 'TOX', 'TIGIT', 'LAG3', 'PDCD1',
               'CTLA4', 'CCL5', 'HOPX', 'RGS1', 'CREM', 'NKG7']

# ══════════════════════════════════════════════════════════════════════
#  Step 1: 加载 H5 → AnnData
# ══════════════════════════════════════════════════════════════════════
print('=' * 60)
print('Step 1: Loading H5 file ...')

H5_FILE = 'data/GSE178341_crc10x_full_c295v4_submit.h5'

with h5py.File(H5_FILE, 'r') as f:
    barcodes = f['matrix/barcodes'][:].astype(str)
    gene_names = f['matrix/features/name'][:].astype(str)
    gene_ids   = f['matrix/features/id'][:].astype(str)
    data    = f['matrix/data'][:]
    indices = f['matrix/indices'][:]
    indptr  = f['matrix/indptr'][:]
    shape   = tuple(f['matrix/shape'][:])  # (n_genes, n_cells)

print(f'  Matrix shape (genes × cells): {shape}')
print(f'  Barcodes: {len(barcodes):,}  Genes: {len(gene_names):,}')

# shape = (n_genes, n_cells), indptr length = n_cells+1 → CSC matrix
# transpose to CSR (cells × genes)
X = csc_matrix((data, indices, indptr), shape=shape).T.tocsr().astype('float32')
print(f'  Sparse matrix: {X.shape}, nnz={X.nnz:,}')

adata = sc.AnnData(
    X   = X,
    obs = pd.DataFrame(index=barcodes),
    var = pd.DataFrame({'gene_ids': gene_ids}, index=gene_names),
)
adata.var_names_make_unique()
print(f'  AnnData: {adata.n_obs:,} cells × {adata.n_vars:,} genes')

# ══════════════════════════════════════════════════════════════════════
#  Step 2: 合并 metatables + cluster 注释
# ══════════════════════════════════════════════════════════════════════
print('\nStep 2: Merging metadata ...')

META_GZ    = 'data/GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz'
CLUSTER_GZ = 'data/GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz'

with gzip.open(META_GZ, 'rt') as f:
    meta = pd.read_csv(f, index_col=0)   # index = cellID
meta.index = meta.index.astype(str)

with gzip.open(CLUSTER_GZ, 'rt') as f:
    cluster = pd.read_csv(f, index_col=0)  # index = sampleID (= cellID)
cluster.index = cluster.index.astype(str)

# 对齐到 adata.obs_names
common_meta    = adata.obs_names.intersection(meta.index)
common_cluster = adata.obs_names.intersection(cluster.index)
print(f'  Meta match:    {len(common_meta):,}/{adata.n_obs:,}')
print(f'  Cluster match: {len(common_cluster):,}/{adata.n_obs:,}')

# 把关键列迁移到 adata.obs
adata.obs = adata.obs.join(
    meta[['SPECIMEN_TYPE', 'MMRStatus', 'TumorStage',
          'HistologicTypeSimple', 'TissueSiteSimple',
          'PID', 'Sex', 'Age', 'NodeStatusSimple',
          'MetastasisStatus', 'MMR_IHC']],
    how='left'
)
adata.obs = adata.obs.join(
    cluster[['clTopLevel', 'clMidwayPr',
             'cl295v11SubShort', 'cl295v11SubFull', 'batchID']],
    how='left'
)

print('  MMRStatus distribution:')
print('   ', adata.obs['MMRStatus'].value_counts().to_dict())
print('  SPECIMEN_TYPE distribution:')
print('   ', adata.obs['SPECIMEN_TYPE'].value_counts().to_dict())

# ══════════════════════════════════════════════════════════════════════
#  Step 3: 筛选 MMRp + Tumor（T）
# ══════════════════════════════════════════════════════════════════════
print('\nStep 3: Filtering to MMRp Tumor cells ...')

mask_mmrp  = adata.obs['MMRStatus'] == 'MMRp'
mask_tumor = adata.obs['SPECIMEN_TYPE'] == 'T'
mask = mask_mmrp & mask_tumor
adata = adata[mask].copy()
print(f'  After filter: {adata.n_obs:,} cells')
print(f'  Patients: {adata.obs["PID"].nunique()}')
print(f'  Cell types (clTopLevel):\n  {adata.obs["clTopLevel"].value_counts().to_dict()}')

# ══════════════════════════════════════════════════════════════════════
#  Step 4: QC filtering
# ══════════════════════════════════════════════════════════════════════
print('\nStep 4: QC filtering ...')

# 标注线粒体基因
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                            log1p=False, inplace=True)

print('  QC stats before filtering:')
print(f'    n_genes_by_counts: median={adata.obs["n_genes_by_counts"].median():.0f}, '
      f'p5={np.percentile(adata.obs["n_genes_by_counts"], 5):.0f}, '
      f'p95={np.percentile(adata.obs["n_genes_by_counts"], 95):.0f}')
print(f'    pct_counts_mt: median={adata.obs["pct_counts_mt"].median():.1f}%, '
      f'p95={np.percentile(adata.obs["pct_counts_mt"], 95):.1f}%')

# 过滤低质量细胞
n_before = adata.n_obs
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)

# 线粒体 < 25%（数据已经是 processed，线粒体比例通常低）
mt_cutoff = min(25, np.percentile(adata.obs['pct_counts_mt'], 97))
adata = adata[adata.obs['pct_counts_mt'] < mt_cutoff].copy()
print(f'  Filtered: {n_before:,} → {adata.n_obs:,} cells ({n_before - adata.n_obs:,} removed)')
print(f'  Genes retained: {adata.n_vars:,}')

# ══════════════════════════════════════════════════════════════════════
#  Step 5: 标准化 + log1p，保留 raw counts
# ══════════════════════════════════════════════════════════════════════
print('\nStep 5: Normalizing ...')
# 存储 raw counts (整型近似，检查数据是否已是浮点处理过的)
vals = adata.X.data[:100]
print(f'  X sample values: {vals[:10]}')
is_counts = np.allclose(vals, np.round(vals), atol=0.01)
print(f'  Looks like raw counts: {is_counts}')

adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ══════════════════════════════════════════════════════════════════════
#  Step 6: HVG + PCA + neighbors + UMAP（快速预览用）
# ══════════════════════════════════════════════════════════════════════
print('\nStep 6: HVG selection + PCA + UMAP ...')
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor='seurat_v3',
                              layer='counts', subset=False)
print(f'  HVGs: {adata.var["highly_variable"].sum():,}')

sc.pp.pca(adata, n_comps=50, use_highly_variable=True)
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=40)
sc.tl.umap(adata, min_dist=0.3, random_state=42)
print('  UMAP done.')

# ══════════════════════════════════════════════════════════════════════
#  Step 7: Pre-driver 基因模块打分
# ══════════════════════════════════════════════════════════════════════
print('\nStep 7: Scoring Pre-driver modules ...')

avail = set(adata.var_names)
tam_use = [g for g in TAM_MODULE if g in avail]
tex_use = [g for g in TEX_MODULE  if g in avail]
print(f'  TAM genes: {len(tam_use)}/{len(TAM_MODULE)} → {tam_use}')
print(f'  Tex genes: {len(tex_use)}/{len(TEX_MODULE)} → {tex_use}')

sc.tl.score_genes(adata, tam_use, score_name='TAM_score',    use_raw=False)
sc.tl.score_genes(adata, tex_use, score_name='CD8Tex_score', use_raw=False)

from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()
tam_norm = scaler.fit_transform(adata.obs[['TAM_score']]).ravel()
tex_norm = scaler.fit_transform(adata.obs[['CD8Tex_score']]).ravel()
adata.obs['Predriver_score'] = np.maximum(tam_norm, tex_norm)

print('\n  Mean scores by cell type (clTopLevel):')
score_df = adata.obs.groupby('clTopLevel')[
    ['TAM_score', 'CD8Tex_score', 'Predriver_score']
].mean().round(3)
print(score_df.to_string())

# ══════════════════════════════════════════════════════════════════════
#  Step 8: 保存
# ══════════════════════════════════════════════════════════════════════
print(f'\nStep 8: Saving to {OUT_H5AD} ...')
adata.write_h5ad(OUT_H5AD, compression='gzip')
sz = os.path.getsize(OUT_H5AD) / 1024**3
print(f'  Saved: {sz:.2f} GB')

# ── 简要统计摘要 ────────────────────────────────────────────────────
print('\n' + '=' * 60)
print('SUMMARY')
print(f'  Total cells (MMRp Tumor): {adata.n_obs:,}')
print(f'  Genes:    {adata.n_vars:,}')
print(f'  Patients: {adata.obs["PID"].nunique()}')
print('\n  Cell type breakdown (clTopLevel):')
print(adata.obs['clTopLevel'].value_counts().to_string())
print('\n  Immune subtypes (clMidwayPr):')
immune_mask = adata.obs['clTopLevel'].isin(['TNKILC', 'Myeloid', 'B'])
print(adata.obs.loc[immune_mask, 'clMidwayPr'].value_counts().to_string())
print('\n  TumorStage:')
print(adata.obs['TumorStage'].value_counts().to_string())
print('\n  Predriver_score mean by clMidwayPr (immune, top 10):')
imm = adata.obs[immune_mask]
print(imm.groupby('clMidwayPr')['Predriver_score'].mean()
       .sort_values(ascending=False).head(10).round(3).to_string())
print('\nDone.')
