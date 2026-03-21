#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GSE144735 KUL3 CRC 数据预处理
- 加载 genes×cells 原始 UMI 矩阵 → AnnData
- 标准化 + log1p
- 用 Pre-driver 基因模块（TAM / CD8_Tex）对每个细胞打分
- 保存 processed_data/GSE144735.h5ad
"""

import scanpy as sc
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import gzip, io, sys, os, warnings
warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

OUT = 'processed_data/GSE144735.h5ad'

# ── 如果已存在则跳过 ──────────────────────────────────────────────────
if os.path.exists(OUT):
    print(f'{OUT} already exists, skipping conversion.')
    adata = sc.read_h5ad(OUT)
else:
    # ── 读取矩阵（genes × cells，tab分隔）────────────────────────────────
    print('Reading expression matrix (genes × cells)...')
    MATRIX_GZ = 'data/GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz'

    # 逐块读取以节省内存
    chunks = []
    gene_names = []
    with gzip.open(MATRIX_GZ, 'rt') as f:
        header = f.readline().strip().split('\t')
        # 第一个字段是行名占位符（空字符串或 'Index'），cell names 从第 1 列开始
        cell_names = header[1:]
        for i, line in enumerate(f):
            parts = line.rstrip('\n').split('\t')
            gene_names.append(parts[0])
            chunks.append([int(x) for x in parts[1:]])
            if (i + 1) % 5000 == 0:
                print(f'  Read {i+1} genes...')

    print(f'  Total genes: {len(gene_names)}, cells: {len(cell_names)}')

    # 构建稀疏矩阵 (cells × genes)
    print('Building sparse matrix...')
    X = csr_matrix(np.array(chunks, dtype=np.int32).T)

    # AnnData: obs=cells, var=genes
    adata = sc.AnnData(X=X,
                       obs=pd.DataFrame(index=cell_names),
                       var=pd.DataFrame(index=gene_names))

    # ── 加载元数据 ────────────────────────────────────────────────────────
    print('Loading annotation...')
    ANN_GZ = 'data/GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz'
    with gzip.open(ANN_GZ, 'rt') as f:
        ann = pd.read_csv(f, sep='\t', index_col=0)

    # 只保留 annotation 中有记录的细胞
    common = adata.obs_names.intersection(ann.index)
    adata = adata[common].copy()
    adata.obs = ann.loc[common]
    print(f'  After annotation filter: {adata.n_obs:,} cells × {adata.n_vars:,} genes')

    # ── 标准化 ────────────────────────────────────────────────────────────
    print('Normalizing...')
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # 保存原始基因符号
    adata.var_names_make_unique()

    # ── 保存 ─────────────────────────────────────────────────────────────
    print(f'Saving to {OUT}...')
    adata.write_h5ad(OUT, compression='gzip')
    print('  Done.')

# ══════════════════════════════════════════════════════════════════════
#  Pre-driver 基因模块打分
# ══════════════════════════════════════════════════════════════════════
print('\nScoring Pre-driver gene modules...')

# TAM 模块（FTH1/SOD2/SOCS3/OSM/HMOX1 + HLA 抗原呈递）
TAM_MODULE = ['FTH1', 'SOD2', 'SOCS3', 'OSM', 'HMOX1', 'FCGR1A',
              'LGALS3', 'CD68', 'HLA-DPB1', 'HLA-DPA1', 'RGS1',
              'CTSB', 'HLA-DRA', 'HLA-DRB1', 'GPR183', 'APOE', 'SPP1']

# CD8_Tex 模块
TEX_MODULE  = ['HAVCR2', 'ENTPD1', 'TOX', 'TIGIT', 'LAG3', 'PDCD1',
               'CTLA4', 'CCL5', 'HOPX', 'RGS1', 'CREM', 'NKG7']

# 仅保留在 var_names 中存在的基因
avail = set(adata.var_names)
tam_use = [g for g in TAM_MODULE if g in avail]
tex_use = [g for g in TEX_MODULE  if g in avail]
print(f'  TAM module genes available: {len(tam_use)}/{len(TAM_MODULE)}  → {tam_use}')
print(f'  Tex module genes available: {len(tex_use)}/{len(TEX_MODULE)}  → {tex_use}')

sc.tl.score_genes(adata, tam_use, score_name='TAM_score',   use_raw=False)
sc.tl.score_genes(adata, tex_use, score_name='CD8Tex_score', use_raw=False)

# Pre-driver 综合评分 = max(TAM, Tex) 标准化到 [0,1]
from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()
tam_norm = scaler.fit_transform(adata.obs[['TAM_score']]).ravel()
tex_norm = scaler.fit_transform(adata.obs[['CD8Tex_score']]).ravel()
adata.obs['Predriver_score'] = np.maximum(tam_norm, tex_norm)

print('\nClass × Cell_type distribution:')
print(pd.crosstab(adata.obs['Class'], adata.obs['Cell_type']))

print('\nMean TAM_score by Class:')
print(adata.obs.groupby('Class')['TAM_score'].mean().round(4))
print('\nMean CD8Tex_score by Class:')
print(adata.obs.groupby('Class')['CD8Tex_score'].mean().round(4))

# 保存打分结果（追加到 h5ad）
print(f'\nUpdating {OUT} with scores...')
adata.write_h5ad(OUT, compression='gzip')
print('Done.')
