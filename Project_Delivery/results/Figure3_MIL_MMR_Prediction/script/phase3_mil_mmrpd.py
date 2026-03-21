#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 3 — 基于注意力多示例学习（ABMIL）的 MMRp/d 异质性建模
GSE178341: 62 patients (34 MMRd + 28 MMRp), T specimens only

Pipeline:
1. 联合加载 MMRp + MMRd 肿瘤细胞 → 共同 PCA 特征空间
2. 构建 Patient-level Bags（每患者 = 一个 bag）
3. 训练 Gated-ABMIL（5-fold CV, cross-validated ROC-AUC）
4. 用最佳模型提取全量细胞注意力分值
5. 统计分析：不同亚群（clMidwayPr）的平均注意力 in MMRp bags
6. 可视化：UMAP 着色注意力 + 亚群注意力小提琴图
"""

import scanpy as sc
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
import h5py, gzip, io, sys, os, warnings, pickle
warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.optim import Adam
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

torch.manual_seed(42)
np.random.seed(42)

plt.rcParams['font.family']  = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42

FIG_DIR = 'figures/phase3'
RES_DIR = 'phase3_results/mil'
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RES_DIR, exist_ok=True)

# ── 超参数 ─────────────────────────────────────────────────────────
N_PCA       = 50       # PCA 维度作为 cell feature
MAX_CELLS   = 2000     # 每 bag 每 epoch 最大采样细胞数（stochastic）
HIDDEN_DIM  = 256      # encoder 隐层宽度
ATT_DIM     = 128      # attention 网络维度
DROPOUT     = 0.25
LR          = 1e-4
WEIGHT_DECAY= 1e-5
EPOCHS      = 80
N_FOLDS     = 5
PATIENCE    = 15       # early stopping patience

# ══════════════════════════════════════════════════════════════════════
#  Step 1: 联合加载 MMRp + MMRd 肿瘤细胞
# ══════════════════════════════════════════════════════════════════════
JOINT_H5AD = 'processed_data/GSE178341_joint.h5ad'

if os.path.exists(JOINT_H5AD):
    print('Loading cached joint data ...')
    adata = sc.read_h5ad(JOINT_H5AD)
    print(f'  {adata.n_obs:,} cells × {adata.n_vars:,} genes')
else:
    print('=' * 60)
    print('Step 1: Loading H5 (T specimens: MMRp + MMRd) ...')

    H5_FILE    = 'data/GSE178341_crc10x_full_c295v4_submit.h5'
    META_GZ    = 'data/GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz'
    CLUSTER_GZ = 'data/GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz'

    with h5py.File(H5_FILE, 'r') as f:
        barcodes   = f['matrix/barcodes'][:].astype(str)
        gene_names = f['matrix/features/name'][:].astype(str)
        gene_ids   = f['matrix/features/id'][:].astype(str)
        data       = f['matrix/data'][:]
        indices    = f['matrix/indices'][:]
        indptr     = f['matrix/indptr'][:]
        shape      = tuple(f['matrix/shape'][:])

    print(f'  Full matrix: {shape}')
    X_full = csc_matrix((data, indices, indptr), shape=shape).T.tocsr().astype('float32')

    with gzip.open(META_GZ, 'rt') as f:
        meta = pd.read_csv(f, index_col=0)
    with gzip.open(CLUSTER_GZ, 'rt') as f:
        cluster = pd.read_csv(f, index_col=0)
    meta.index    = meta.index.astype(str)
    cluster.index = cluster.index.astype(str)

    # 构建临时 AnnData
    adata = sc.AnnData(
        X   = X_full,
        obs = pd.DataFrame(index=barcodes),
        var = pd.DataFrame({'gene_ids': gene_ids}, index=gene_names),
    )
    adata.var_names_make_unique()

    # 合并 meta
    adata.obs = adata.obs.join(
        meta[['SPECIMEN_TYPE', 'MMRStatus', 'TumorStage',
              'PID', 'Sex', 'Age', 'NodeStatusSimple']],
        how='left'
    )
    adata.obs = adata.obs.join(
        cluster[['clTopLevel', 'clMidwayPr', 'cl295v11SubShort',
                 'cl295v11SubFull', 'batchID']],
        how='left'
    )

    # 仅保留 T 标本 + MMRp/MMRd
    mask = (adata.obs['SPECIMEN_TYPE'] == 'T') & \
           (adata.obs['MMRStatus'].isin(['MMRp', 'MMRd']))
    adata = adata[mask].copy()
    print(f'  After T+MMR filter: {adata.n_obs:,} cells')
    print(f'  MMRp: {(adata.obs["MMRStatus"]=="MMRp").sum():,}  '
          f'MMRd: {(adata.obs["MMRStatus"]=="MMRd").sum():,}')

    # QC
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                               log1p=False, inplace=True)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=10)
    mt_cutoff = min(25, np.percentile(adata.obs['pct_counts_mt'], 97))
    adata = adata[adata.obs['pct_counts_mt'] < mt_cutoff].copy()
    print(f'  After QC: {adata.n_obs:,} cells, {adata.n_vars:,} genes')

    # 标准化
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # HVG + PCA（共同特征空间）
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor='seurat_v3',
                                 layer='counts', subset=False)
    sc.pp.pca(adata, n_comps=N_PCA, use_highly_variable=True)

    # UMAP（仅供可视化）
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=40)
    sc.tl.umap(adata, min_dist=0.3, random_state=42)

    adata.write_h5ad(JOINT_H5AD, compression='gzip')
    print(f'  Saved: {JOINT_H5AD}')

print(f'\nMMRp patients: {adata.obs[adata.obs["MMRStatus"]=="MMRp"]["PID"].nunique()}')
print(f'MMRd patients: {adata.obs[adata.obs["MMRStatus"]=="MMRd"]["PID"].nunique()}')
print('Cell type distribution:')
print(adata.obs.groupby('MMRStatus')['clTopLevel'].value_counts().unstack(fill_value=0))

# ══════════════════════════════════════════════════════════════════════
#  Step 2: 构建 Bag Dataset
# ══════════════════════════════════════════════════════════════════════
print('\nStep 2: Building bags ...')

pca_feat = adata.obsm['X_pca']          # (N_cells, N_PCA)
mmr_label = (adata.obs['MMRStatus'] == 'MMRd').astype(int).values  # 0=MMRp, 1=MMRd
pids = adata.obs['PID'].values

# patient → (cell indices, label)
unique_pids  = adata.obs['PID'].unique()
patient_bags = {}
for pid in unique_pids:
    idx = np.where(pids == pid)[0]
    lbl = mmr_label[idx[0]]
    patient_bags[pid] = {'idx': idx, 'label': lbl}

all_pids  = np.array(list(patient_bags.keys()))
all_labels = np.array([patient_bags[p]['label'] for p in all_pids])
print(f'  Bags: {len(all_pids)} patients '
      f'(MMRd={all_labels.sum()}, MMRp={(1-all_labels).sum()})')
print(f'  Cells per bag: median={np.median([len(patient_bags[p]["idx"]) for p in all_pids]):.0f}')

# ══════════════════════════════════════════════════════════════════════
#  Step 3: ABMIL 模型
# ══════════════════════════════════════════════════════════════════════
class GatedABMIL(nn.Module):
    """Gated Attention-Based MIL (Ilse et al. 2018)"""
    def __init__(self, input_dim=N_PCA, hidden_dim=HIDDEN_DIM,
                 att_dim=ATT_DIM, dropout=DROPOUT):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
        )
        # Gated attention
        self.att_V = nn.Sequential(nn.Linear(hidden_dim, att_dim), nn.Tanh())
        self.att_U = nn.Sequential(nn.Linear(hidden_dim, att_dim), nn.Sigmoid())
        self.att_w = nn.Linear(att_dim, 1, bias=False)
        # Bag classifier
        self.classifier = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 1),
        )

    def forward(self, x, return_attention=False):
        """x: (N_cells, input_dim)"""
        H = self.encoder(x)                         # (N, hidden)
        A = self.att_w(self.att_V(H) * self.att_U(H))  # (N, 1)
        A = F.softmax(A, dim=0)                     # normalize over cells
        M = (A * H).sum(dim=0, keepdim=True)        # (1, hidden) bag repr
        logit = self.classifier(M).squeeze()         # scalar
        if return_attention:
            return logit, A.squeeze()                # also return per-cell A
        return logit


def train_epoch(model, optimizer, bags_train, pca_feat, patient_bags,
                max_cells=MAX_CELLS):
    model.train()
    total_loss = 0
    pids_shuf = bags_train.copy()
    np.random.shuffle(pids_shuf)
    for pid in pids_shuf:
        bag_info = patient_bags[pid]
        idx = bag_info['idx']
        label = torch.tensor(bag_info['label'], dtype=torch.float32)
        # stochastic bag subsampling
        if len(idx) > max_cells:
            idx = np.random.choice(idx, max_cells, replace=False)
        x = torch.tensor(pca_feat[idx], dtype=torch.float32)
        optimizer.zero_grad()
        logit = model(x)
        loss = F.binary_cross_entropy_with_logits(logit, label)
        loss.backward()
        optimizer.step()
        total_loss += loss.item()
    return total_loss / len(pids_shuf)


@torch.no_grad()
def eval_bags(model, pids_eval, pca_feat, patient_bags, max_cells=None):
    model.eval()
    probs, labels = [], []
    for pid in pids_eval:
        bag_info = patient_bags[pid]
        idx = bag_info['idx']
        if max_cells and len(idx) > max_cells:
            # use all cells for eval
            pass
        x = torch.tensor(pca_feat[idx], dtype=torch.float32)
        logit = model(x)
        probs.append(torch.sigmoid(logit).item())
        labels.append(bag_info['label'])
    return np.array(probs), np.array(labels)


# ══════════════════════════════════════════════════════════════════════
#  Step 4: 5-Fold CV 训练
# ══════════════════════════════════════════════════════════════════════
print('\nStep 4: 5-Fold Cross-Validation ...')

skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=42)
fold_aucs  = []
best_models = {}

for fold, (train_idx, val_idx) in enumerate(skf.split(all_pids, all_labels)):
    pids_train = all_pids[train_idx]
    pids_val   = all_pids[val_idx]
    print(f'  Fold {fold+1}/{N_FOLDS}: train={len(pids_train)}, val={len(pids_val)}')

    model = GatedABMIL()
    opt   = Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)

    best_auc   = 0
    best_state = None
    patience_cnt = 0

    for epoch in range(EPOCHS):
        loss = train_epoch(model, opt, list(pids_train), pca_feat, patient_bags)
        probs_val, labels_val = eval_bags(model, pids_val, pca_feat, patient_bags)
        if len(np.unique(labels_val)) > 1:
            auc = roc_auc_score(labels_val, probs_val)
        else:
            auc = 0.5
        if auc > best_auc:
            best_auc   = auc
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
            patience_cnt = 0
        else:
            patience_cnt += 1
        if patience_cnt >= PATIENCE:
            print(f'    Early stop at epoch {epoch+1}, best AUC={best_auc:.4f}')
            break

    if epoch == EPOCHS - 1:
        print(f'    Completed {EPOCHS} epochs, best AUC={best_auc:.4f}')

    fold_aucs.append(best_auc)
    best_models[fold] = best_state

print(f'\n  CV AUC: {[f"{a:.4f}" for a in fold_aucs]}')
print(f'  Mean AUC: {np.mean(fold_aucs):.4f} ± {np.std(fold_aucs):.4f}')

# 保存 CV 结果
cv_result = pd.DataFrame({'fold': range(1, N_FOLDS+1), 'auc': fold_aucs})
cv_result.to_csv(f'{RES_DIR}/cv_auc.csv', index=False)

# ══════════════════════════════════════════════════════════════════════
#  Step 5: 用全部患者重新训练最终模型（用最好 fold 的超参）
#  提取全量细胞注意力分值
# ══════════════════════════════════════════════════════════════════════
print('\nStep 5: Training final model on all patients & extracting attention ...')

# 选择 AUC 最好的 fold 的 model state 作为初始化
best_fold = int(np.argmax(fold_aucs))
print(f'  Best fold: {best_fold+1} (AUC={fold_aucs[best_fold]:.4f})')

final_model = GatedABMIL()
final_model.load_state_dict(best_models[best_fold])

# Fine-tune 20 epochs on all data
opt_final = Adam(final_model.parameters(), lr=LR * 0.5, weight_decay=WEIGHT_DECAY)
for epoch in range(20):
    loss = train_epoch(final_model, opt_final, list(all_pids),
                       pca_feat, patient_bags)

torch.save(final_model.state_dict(), f'{RES_DIR}/abmil_final.pt')

# ── 提取每个细胞的注意力分值 ────────────────────────────────────────
print('  Extracting per-cell attention scores ...')
final_model.eval()
att_scores = np.zeros(adata.n_obs, dtype=np.float32)

with torch.no_grad():
    for pid in all_pids:
        idx = patient_bags[pid]['idx']
        x   = torch.tensor(pca_feat[idx], dtype=torch.float32)
        _, A = final_model(x, return_attention=True)
        att_scores[idx] = A.numpy()

adata.obs['attention_score'] = att_scores
print(f'  Attention score range: [{att_scores.min():.6f}, {att_scores.max():.6f}]')
print(f'  Attention score mean by MMR:')
print(adata.obs.groupby('MMRStatus')['attention_score'].mean())

# ══════════════════════════════════════════════════════════════════════
#  Step 6: 统计分析 — 各亚群注意力分值
# ══════════════════════════════════════════════════════════════════════
print('\nStep 6: Attention analysis by cell subtype ...')

# 分别分析 MMRp 和 MMRd 患者中各亚群的注意力
def analyze_attention_by_subtype(adata_sub, subtype_col, mmr_status, top_n=15):
    sub = adata_sub[adata_sub.obs['MMRStatus'] == mmr_status]
    grp = sub.obs.groupby(subtype_col)['attention_score']
    stats = grp.agg(['mean', 'median', 'std', 'count']).reset_index()
    stats.columns = [subtype_col, 'mean_att', 'median_att', 'std_att', 'n_cells']
    stats = stats[stats['n_cells'] >= 50].sort_values('mean_att', ascending=False)
    return stats.head(top_n)

# 用 clMidwayPr（作者标注的中等粒度亚群）
stats_mmrp = analyze_attention_by_subtype(adata, 'clMidwayPr', 'MMRp')
stats_mmrd = analyze_attention_by_subtype(adata, 'clMidwayPr', 'MMRd')

print('\nMMRp — top attention subtypes:')
print(stats_mmrp.to_string(index=False))
print('\nMMRd — top attention subtypes:')
print(stats_mmrd.to_string(index=False))

stats_mmrp.to_csv(f'{RES_DIR}/attention_by_subtype_MMRp.csv', index=False)
stats_mmrd.to_csv(f'{RES_DIR}/attention_by_subtype_MMRd.csv', index=False)

# Mann-Whitney U: Macro vs all others (in MMRp bags)
mmrp_cells = adata.obs[adata.obs['MMRStatus'] == 'MMRp']
macro_att = mmrp_cells[mmrp_cells['clMidwayPr'] == 'Macro']['attention_score'].values
tcd8_att  = mmrp_cells[mmrp_cells['clMidwayPr'] == 'TCD8']['attention_score'].values
other_att = mmrp_cells[~mmrp_cells['clMidwayPr'].isin(['Macro','TCD8'])]['attention_score'].values

mw_macro, p_macro = mannwhitneyu(macro_att, other_att, alternative='greater')
mw_tcd8,  p_tcd8  = mannwhitneyu(tcd8_att,  other_att, alternative='greater')
print(f'\nMann-Whitney (MMRp bags):')
print(f'  Macro vs Others: U={mw_macro:.0f}, p={p_macro:.4e}')
print(f'  TCD8  vs Others: U={mw_tcd8:.0f},  p={p_tcd8:.4e}')

# ══════════════════════════════════════════════════════════════════════
#  Step 7: 可视化
# ══════════════════════════════════════════════════════════════════════
print('\nStep 7: Plotting ...')

# ── 颜色映射 ────────────────────────────────────────────────────────
TOP_LEVEL_COLORS = {
    'TNKILC':  '#3C5488', 'Myeloid': '#E64B35', 'B':    '#00A087',
    'Plasma':  '#F39B7F', 'Epi':     '#8491B4', 'Strom':'#7E6148',
    'Mast':    '#B09C85',
}
MMR_COLORS = {'MMRp': '#E64B35', 'MMRd': '#4DBBD5'}

xy = adata.obsm['X_umap']

# ── Figure layout: 3×3 grid ─────────────────────────────────────────
fig = plt.figure(figsize=(20, 18))
fig.patch.set_facecolor('white')
gs = gridspec.GridSpec(3, 3, figure=fig,
                       hspace=0.45, wspace=0.35,
                       left=0.06, right=0.97, top=0.93, bottom=0.06)

# ── 1. UMAP — MMR status ─────────────────────────────────────────────
ax = fig.add_subplot(gs[0, 0])
for mmr, color in MMR_COLORS.items():
    m = adata.obs['MMRStatus'] == mmr
    ax.scatter(xy[m, 0], xy[m, 1], c=color, s=0.8, alpha=0.4,
               linewidths=0, rasterized=True, label=mmr)
ax.legend(fontsize=8, frameon=False, markerscale=4, loc='upper right')
ax.set_title('MMR Status', fontsize=10, fontweight='bold')
ax.set_xticks([]); ax.set_yticks([])
ax.spines[['top','right','bottom','left']].set_visible(False)

# ── 2. UMAP — Cell type ─────────────────────────────────────────────
ax = fig.add_subplot(gs[0, 1])
for ct, color in TOP_LEVEL_COLORS.items():
    m = adata.obs['clTopLevel'] == ct
    ax.scatter(xy[m, 0], xy[m, 1], c=color, s=0.8, alpha=0.4,
               linewidths=0, rasterized=True, label=ct)
ax.legend(fontsize=7, frameon=False, markerscale=3, loc='upper right',
           ncol=1)
ax.set_title('Cell Type (clTopLevel)', fontsize=10, fontweight='bold')
ax.set_xticks([]); ax.set_yticks([])
ax.spines[['top','right','bottom','left']].set_visible(False)

# ── 3. UMAP — Attention score (全体) ────────────────────────────────
ax = fig.add_subplot(gs[0, 2])
vmax_att = np.percentile(att_scores, 98)
sc3 = ax.scatter(xy[:, 0], xy[:, 1], c=att_scores,
                  cmap='YlOrRd', vmin=0, vmax=vmax_att,
                  s=0.8, alpha=0.6, linewidths=0, rasterized=True)
plt.colorbar(sc3, ax=ax, shrink=0.7, pad=0.02, label='Attention score')
ax.set_title('ABMIL Attention Score', fontsize=10, fontweight='bold')
ax.set_xticks([]); ax.set_yticks([])
ax.spines[['top','right','bottom','left']].set_visible(False)

# ── 4. UMAP — Attention score (MMRp cells only) ──────────────────────
ax = fig.add_subplot(gs[1, 0])
m_mmrp = adata.obs['MMRStatus'] == 'MMRp'
m_mmrd = adata.obs['MMRStatus'] == 'MMRd'
# Gray background: MMRd cells
ax.scatter(xy[m_mmrd, 0], xy[m_mmrd, 1], c='#dddddd',
           s=0.6, alpha=0.3, linewidths=0, rasterized=True)
sc4 = ax.scatter(xy[m_mmrp, 0], xy[m_mmrp, 1],
                  c=att_scores[m_mmrp],
                  cmap='YlOrRd', vmin=0, vmax=vmax_att,
                  s=1.0, alpha=0.7, linewidths=0, rasterized=True)
plt.colorbar(sc4, ax=ax, shrink=0.7, pad=0.02, label='Attention (MMRp)')
ax.set_title('Attention Score — MMRp cells', fontsize=10, fontweight='bold')
ax.set_xticks([]); ax.set_yticks([])
ax.spines[['top','right','bottom','left']].set_visible(False)

# ── 5. UMAP — Attention score (MMRd cells only) ──────────────────────
ax = fig.add_subplot(gs[1, 1])
ax.scatter(xy[m_mmrp, 0], xy[m_mmrp, 1], c='#dddddd',
           s=0.6, alpha=0.3, linewidths=0, rasterized=True)
sc5 = ax.scatter(xy[m_mmrd, 0], xy[m_mmrd, 1],
                  c=att_scores[m_mmrd],
                  cmap='YlGnBu', vmin=0, vmax=vmax_att,
                  s=1.0, alpha=0.7, linewidths=0, rasterized=True)
plt.colorbar(sc5, ax=ax, shrink=0.7, pad=0.02, label='Attention (MMRd)')
ax.set_title('Attention Score — MMRd cells', fontsize=10, fontweight='bold')
ax.set_xticks([]); ax.set_yticks([])
ax.spines[['top','right','bottom','left']].set_visible(False)

# ── 6. CV ROC AUC bar ───────────────────────────────────────────────
ax = fig.add_subplot(gs[1, 2])
fold_labels = [f'Fold {i+1}' for i in range(N_FOLDS)]
colors_bar  = ['#E64B35' if a == max(fold_aucs) else '#8491B4'
               for a in fold_aucs]
bars = ax.bar(fold_labels, fold_aucs, color=colors_bar, alpha=0.85,
              edgecolor='white', width=0.6)
ax.axhline(np.mean(fold_aucs), color='black', lw=1.5, ls='--',
           label=f'Mean={np.mean(fold_aucs):.3f}')
ax.axhline(0.5, color='grey', lw=0.8, ls=':', alpha=0.7)
for bar, auc in zip(bars, fold_aucs):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
            f'{auc:.3f}', ha='center', va='bottom', fontsize=8)
ax.set_ylim(0.4, 1.05)
ax.set_ylabel('ROC-AUC', fontsize=9)
ax.set_title(f'5-Fold CV Performance\nMMRp vs MMRd',
              fontsize=10, fontweight='bold')
ax.legend(fontsize=8, frameon=False)
ax.spines[['top', 'right']].set_visible(False)
ax.set_xticklabels(fold_labels, fontsize=8.5)

# ── 7. 亚群注意力 violin（免疫细胞，MMRp bags）──────────────────────
ax = fig.add_subplot(gs[2, :2])

# 筛选免疫细胞（T + NK + Myeloid）
immune_types = ['TNKILC', 'Myeloid']
imm_mask = (adata.obs['MMRStatus'] == 'MMRp') & \
           (adata.obs['clTopLevel'].isin(immune_types))
imm_df = adata.obs[imm_mask][['clMidwayPr', 'attention_score', 'clTopLevel']].copy()

# 按亚群均值排序
subtype_order = (imm_df.groupby('clMidwayPr')['attention_score']
                 .mean().sort_values(ascending=False).index.tolist())
# 过滤太少的亚群
subtype_order = [s for s in subtype_order
                 if (imm_df['clMidwayPr'] == s).sum() >= 50]

positions = range(len(subtype_order))
subtype_colors = {
    'Macro': '#E64B35', 'TCD8': '#3C5488', 'Mono': '#F39B7F',
    'DC': '#E64B35', 'TCD4': '#3C5488', 'NK': '#00A087',
    'B': '#00A087', 'Tgd': '#8491B4', 'TZBTB16': '#8491B4',
    'Mast': '#B09C85', 'Granulo': '#B09C85', 'ILC': '#7E6148',
}

parts = ax.violinplot(
    [imm_df[imm_df['clMidwayPr'] == s]['attention_score'].values
     for s in subtype_order],
    positions=list(positions), widths=0.7,
    showmedians=True, showextrema=False
)
for i, (body, st) in enumerate(zip(parts['bodies'], subtype_order)):
    col = subtype_colors.get(st, '#aaaaaa')
    body.set_facecolor(col)
    body.set_alpha(0.6)
parts['cmedians'].set_color('black')
parts['cmedians'].set_linewidth(1.5)

# 在 Macro 和 TCD8 上方标注 p 值
def sig_label(p):
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return 'ns'

y_max = imm_df['attention_score'].quantile(0.995)
for st, p_val in [('Macro', p_macro), ('TCD8', p_tcd8)]:
    if st in subtype_order:
        xi = subtype_order.index(st)
        ax.text(xi, y_max * 1.02, sig_label(p_val),
                ha='center', fontsize=10, fontweight='bold',
                color=subtype_colors.get(st, 'black'))

ax.set_xticks(list(positions))
ax.set_xticklabels(subtype_order, rotation=35, ha='right', fontsize=9)
ax.set_ylabel('Attention Score', fontsize=9)
ax.set_title('Per-cell Attention Score by Immune Subtype (MMRp bags)',
              fontsize=10, fontweight='bold')
ax.spines[['top', 'right']].set_visible(False)

# 添加 Pre-driver 区域标注
for st in ['Macro', 'TCD8']:
    if st in subtype_order:
        xi = subtype_order.index(st)
        ax.axvspan(xi - 0.4, xi + 0.4, alpha=0.08, color='#B22222', zorder=0)
        ax.text(xi, -0.00008, 'Pre-driver', ha='center', va='top',
                fontsize=7, color='#B22222', style='italic')

# ── 8. 亚群 mean attention heatmap（MMRp vs MMRd 对比）───────────────
ax = fig.add_subplot(gs[2, 2])

immune_types_all = ['TNKILC', 'Myeloid']
imm_df_all = adata.obs[adata.obs['clTopLevel'].isin(immune_types_all)][
    ['clMidwayPr', 'attention_score', 'MMRStatus']
].copy()

pivot = imm_df_all.groupby(['clMidwayPr', 'MMRStatus'])['attention_score'].mean().unstack()
pivot = pivot[['MMRp', 'MMRd']].dropna()
pivot = pivot.loc[pivot.count(axis=1) == 2]
pivot = pivot.sort_values('MMRp', ascending=False).head(12)

im = ax.imshow(pivot.values.T, aspect='auto', cmap='RdBu_r',
               vmin=pivot.values.min(), vmax=pivot.values.max())
plt.colorbar(im, ax=ax, shrink=0.6, pad=0.04, label='Mean attention')
ax.set_yticks([0, 1])
ax.set_yticklabels(['MMRp', 'MMRd'], fontsize=9)
ax.set_xticks(range(len(pivot)))
ax.set_xticklabels(pivot.index, rotation=40, ha='right', fontsize=8.5)
ax.set_title('Mean Attention by Subtype\n(MMRp vs MMRd)', fontsize=9,
              fontweight='bold')
ax.spines[['top','right','bottom','left']].set_visible(False)

# ── Global title ────────────────────────────────────────────────────
fig.text(0.5, 0.975,
         f'Figure 3A — ABMIL Attention Analysis: MMRp vs MMRd  '
         f'(5-fold CV AUC={np.mean(fold_aucs):.3f}±{np.std(fold_aucs):.3f})',
         fontsize=13, fontweight='bold', ha='center', va='top')

out_name = 'figure3A_mil_attention'
fig.savefig(f'{FIG_DIR}/{out_name}.png', dpi=200, bbox_inches='tight', facecolor='white')
fig.savefig(f'{FIG_DIR}/{out_name}.pdf',           bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'  Saved: {FIG_DIR}/{out_name}.png/.pdf')

# ── 保存注意力分值到 adata ─────────────────────────────────────────
print('\nSaving updated adata with attention scores ...')
adata.write_h5ad(JOINT_H5AD, compression='gzip')
print(f'  Updated: {JOINT_H5AD}')

# 保存详细结果表
att_df = adata.obs[['PID', 'MMRStatus', 'clTopLevel', 'clMidwayPr',
                     'cl295v11SubShort', 'attention_score']].copy()
att_df.to_csv(f'{RES_DIR}/cell_attention_scores.csv')
print(f'  Saved: {RES_DIR}/cell_attention_scores.csv')

print('\n' + '=' * 60)
print('DONE')
print(f'  CV AUC: {np.mean(fold_aucs):.4f} ± {np.std(fold_aucs):.4f}')
print(f'  Macro attention (MMRp) vs others: p={p_macro:.4e}')
print(f'  TCD8  attention (MMRp) vs others: p={p_tcd8:.4e}')
