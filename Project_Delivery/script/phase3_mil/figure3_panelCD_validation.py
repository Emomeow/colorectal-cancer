#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 3 Panel C & D — Weight Transfer Validation
训练集：GSE178341 免疫细胞（MMRp/d 标签），以 766 共同基因 PCA 为特征
验证集：GSE132465 免疫细胞（冻结权重推断，提取 Attention Score）

Panel C: Attention Score UMAP（GSE132465，高分区域验证 Pre-driver 亚群定位）
Panel D: 各亚型平均 Attention 权重排名（发现集 vs 验证集对比）
"""

import scanpy as sc
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.optim import Adam
from sklearn.decomposition import PCA
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import normalize
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import io, sys, os, warnings
warnings.filterwarnings('ignore')

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

torch.manual_seed(42)
np.random.seed(42)

plt.rcParams.update({
    'font.family':  'Arial',
    'pdf.fonttype': 42,
    'ps.fonttype':  42,
})

FIG_DIR = 'figures/phase3'
RES_DIR = 'phase3_results/transfer'
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RES_DIR, exist_ok=True)

# ── 超参数 ─────────────────────────────────────────────────────────────────
N_COMPS   = 50
MAX_CELLS = 1500
HIDDEN_DIM = 256
ATT_DIM    = 128
DROPOUT    = 0.25
LR         = 1e-4
WD         = 1e-5
EPOCHS     = 100
N_FOLDS    = 5
PATIENCE   = 20

PRE_DRIVER = {'Macro', 'TCD8'}
COLORS_MMR = {'MMRp': '#E64B35', 'MMRd': '#4DBBD5'}

# GSE132465 亚型 → 标准化名称映射
GSE132_SUBTYPE_MAP = {
    'CD8+ T cells':              'TCD8',
    'CD4+ T cells':              'TCD4',
    'T helper 17 cells':         'TCD4',
    'T follicular helper cells': 'TCD4',
    'Regulatory T cells':        'Treg',
    'gamma delta T cells':       'Tgd',
    'NK cells':                  'NK',
    'SPP1+':                     'Macro',
    'Pro-inflammatory':          'Macro',
    'cDC':                       'DC',
    'Proliferating':             'Other',
    'Unknown':                   'Other',
}
SUBTYPE_COLORS = {
    'Macro': '#E64B35', 'TCD8': '#3C5488', 'TCD4': '#00A087',
    'Treg': '#F39B7F', 'NK': '#8491B4', 'DC': '#91D1C2',
    'Tgd': '#B09C85', 'Other': '#DDDDDD',
}

# ══════════════════════════════════════════════════════════════════════════
#  Step 1: 加载 GSE178341 免疫数据，确定共同 HVG
# ══════════════════════════════════════════════════════════════════════════
print('=' * 65)
print('Step 1: Loading GSE178341 immune data ...')
adata178 = sc.read_h5ad('processed_data/GSE178341_immune_mil.h5ad')
print(f'  {adata178.n_obs:,} cells × {adata178.n_vars:,} genes')

# 加载 merged 数据获取基因列表
print('  Loading merged gene list ...')
merged_genes = set(sc.read_h5ad('processed_data/merged_annotated.h5ad', backed='r').var_names)

# 共同 HVG
hvg_all = adata178.var_names[adata178.var['highly_variable']].tolist()
common_hvg = [g for g in hvg_all if g in merged_genes]
print(f'  HVGs from GSE178341: {len(hvg_all)}, common with merged: {len(common_hvg)}')

# ══════════════════════════════════════════════════════════════════════════
#  Step 2: 在 GSE178341 上拟合 PCA（common HVG 特征空间）
# ══════════════════════════════════════════════════════════════════════════
print('\nStep 2: Fitting PCA on GSE178341 common HVGs ...')

import scipy.sparse as sp
X178 = adata178[:, common_hvg].X
if sp.issparse(X178):
    X178 = X178.toarray()
X178 = X178.astype(np.float32)
print(f'  X178 shape: {X178.shape}')

pca_model = PCA(n_components=N_COMPS, svd_solver='arpack', random_state=42)
X178_pca = pca_model.fit_transform(X178)
print(f'  PCA fitted: explained variance ratio sum = {pca_model.explained_variance_ratio_.sum():.3f}')

# Bag 构建
pids178     = adata178.obs['PID'].values
mmr178      = (adata178.obs['MMRStatus'] == 'MMRd').astype(int).values
unique_pids = adata178.obs['PID'].unique()
bags178 = {}
for pid in unique_pids:
    idx = np.where(pids178 == pid)[0]
    bags178[pid] = {'idx': idx, 'label': int(mmr178[idx[0]])}

all_pids178   = np.array(list(bags178.keys()))
all_labels178 = np.array([bags178[p]['label'] for p in all_pids178])
print(f'  Bags: {len(all_pids178)} (MMRd={all_labels178.sum()}, MMRp={(1-all_labels178).sum()})')

# ══════════════════════════════════════════════════════════════════════════
#  Step 3: ABMIL 模型定义 + 训练
# ══════════════════════════════════════════════════════════════════════════
class GatedABMIL(nn.Module):
    def __init__(self):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(N_COMPS, HIDDEN_DIM), nn.ReLU(), nn.Dropout(DROPOUT),
            nn.Linear(HIDDEN_DIM, HIDDEN_DIM), nn.ReLU(),
        )
        self.att_V = nn.Sequential(nn.Linear(HIDDEN_DIM, ATT_DIM), nn.Tanh())
        self.att_U = nn.Sequential(nn.Linear(HIDDEN_DIM, ATT_DIM), nn.Sigmoid())
        self.att_w = nn.Linear(ATT_DIM, 1, bias=False)
        self.classifier = nn.Sequential(nn.Dropout(DROPOUT), nn.Linear(HIDDEN_DIM, 1))

    def forward(self, x, return_attention=False):
        H = self.encoder(x)
        A = F.softmax(self.att_w(self.att_V(H) * self.att_U(H)), dim=0)
        M = (A * H).sum(dim=0, keepdim=True)
        logit = self.classifier(M).squeeze()
        if return_attention:
            return logit, A.squeeze()
        return logit


def train_epoch(model, opt, pids_list, X_pca, bags):
    model.train()
    order = pids_list.copy()
    np.random.shuffle(order)
    for pid in order:
        idx = bags[pid]['idx']
        if len(idx) > MAX_CELLS:
            idx = np.random.choice(idx, MAX_CELLS, replace=False)
        x   = torch.tensor(X_pca[idx], dtype=torch.float32)
        lbl = torch.tensor(bags[pid]['label'], dtype=torch.float32)
        opt.zero_grad()
        F.binary_cross_entropy_with_logits(model(x), lbl).backward()
        opt.step()


@torch.no_grad()
def predict(model, pids_list, X_pca, bags):
    model.eval()
    probs, labels = [], []
    for pid in pids_list:
        idx = bags[pid]['idx']
        probs.append(torch.sigmoid(model(torch.tensor(X_pca[idx], dtype=torch.float32))).item())
        labels.append(bags[pid]['label'])
    return np.array(probs), np.array(labels)


print('\nStep 3: Training ABMIL on common HVG PCA ...')
skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=42)
fold_aucs, best_states = [], {}

for fold, (tr_idx, va_idx) in enumerate(skf.split(all_pids178, all_labels178)):
    pids_tr = all_pids178[tr_idx].tolist()
    pids_va = all_pids178[va_idx].tolist()
    model = GatedABMIL()
    opt   = Adam(model.parameters(), lr=LR, weight_decay=WD)
    best_auc, best_state, patience_cnt = 0.0, None, 0
    for epoch in range(EPOCHS):
        train_epoch(model, opt, pids_tr, X178_pca, bags178)
        pv, lv = predict(model, pids_va, X178_pca, bags178)
        auc = roc_auc_score(lv, pv) if len(np.unique(lv)) > 1 else 0.5
        if auc > best_auc:
            best_auc, best_state, patience_cnt = auc, {k: v.clone() for k, v in model.state_dict().items()}, 0
        else:
            patience_cnt += 1
        if patience_cnt >= PATIENCE:
            break
    fold_aucs.append(best_auc)
    best_states[fold] = best_state
    print(f'  Fold {fold+1}: AUC={best_auc:.4f}')

print(f'  Mean CV AUC: {np.mean(fold_aucs):.4f} ± {np.std(fold_aucs):.4f}')

# 最终模型（在全量 GSE178341 上微调）
best_fold = int(np.argmax(fold_aucs))
final_model = GatedABMIL()
final_model.load_state_dict(best_states[best_fold])
opt_f = Adam(final_model.parameters(), lr=LR * 0.5, weight_decay=WD)
for _ in range(30):
    train_epoch(final_model, opt_f, list(all_pids178), X178_pca, bags178)
torch.save(final_model.state_dict(), f'{RES_DIR}/abmil_transfer.pt')
print('  Final model saved.')

# ── 提取 GSE178341 per-cell attention（MMRp 患者，用于 Panel D 对比）────────
print('  Extracting GSE178341 attention scores ...')
final_model.eval()
att178 = np.zeros(adata178.n_obs, dtype=np.float32)
with torch.no_grad():
    for pid in all_pids178:
        idx = bags178[pid]['idx']
        x   = torch.tensor(X178_pca[idx], dtype=torch.float32)
        _, A = final_model(x, return_attention=True)
        att178[idx] = A.numpy()

adata178.obs['attention_transfer'] = att178
# 按亚型统计（MMRp 患者，与 panel D 对比）
mmrp178_obs = adata178.obs[adata178.obs['MMRStatus'] == 'MMRp']
att178_by_ct = mmrp178_obs.groupby('clMidwayPr')['attention_transfer'].agg(
    ['mean', 'median', 'std', 'count']).reset_index()
att178_by_ct.columns = ['subtype', 'mean', 'median', 'std', 'n']
att178_by_ct = att178_by_ct[att178_by_ct['n'] >= 50].sort_values('mean', ascending=False)
print(f'  GSE178341 attention by subtype (MMRp):\n{att178_by_ct.to_string(index=False)}')

# ══════════════════════════════════════════════════════════════════════════
#  Step 4: 加载 GSE132465 免疫细胞，投影到共同 PCA 空间
# ══════════════════════════════════════════════════════════════════════════
print('\nStep 4: Processing GSE132465 immune cells ...')
adata_merged = sc.read_h5ad('processed_data/merged_annotated.h5ad')
# 筛选 GSE132465 肿瘤免疫细胞
IMMUNE_TYPES132 = ['T cells', 'Myeloids', 'NK cells']
mask_gse = ((adata_merged.obs['batch'] == 'GSE132465') &
             adata_merged.obs['unified_cell_type'].isin(IMMUNE_TYPES132) &
             (adata_merged.obs['Class'] == 'Tumor'))
adata132 = adata_merged[mask_gse].copy()
del adata_merged  # 释放内存
print(f'  GSE132465 tumor immune cells: {adata132.n_obs:,}')
print(f'  Patients: {adata132.obs["Patient"].nunique()}')

# 标准化（GSE132465 .X 为原始 counts）
print('  Normalizing GSE132465 ...')
sc.pp.normalize_total(adata132, target_sum=1e4)
sc.pp.log1p(adata132)

# 构建 common HVG 特征矩阵（缺失基因用 0）
print('  Building feature matrix ...')
X132 = np.zeros((adata132.n_obs, len(common_hvg)), dtype=np.float32)
gse132_gene_idx = {g: i for i, g in enumerate(adata132.var_names)}
for j, gene in enumerate(common_hvg):
    if gene in gse132_gene_idx:
        col_idx = gse132_gene_idx[gene]
        col_data = adata132.X[:, col_idx]
        if sp.issparse(col_data):
            col_data = col_data.toarray().ravel()
        else:
            col_data = np.asarray(col_data).ravel()
        X132[:, j] = col_data

# 投影到 GSE178341 PCA 空间
X132_pca = pca_model.transform(X132).astype(np.float32)
print(f'  X132_pca shape: {X132_pca.shape}')

# UMAP for GSE132465
print('  Computing UMAP for GSE132465 ...')
from sklearn.neighbors import NearestNeighbors
import umap as umap_lib
reducer = umap_lib.UMAP(n_components=2, n_neighbors=20, min_dist=0.3, random_state=42)
X132_umap = reducer.fit_transform(X132_pca)
adata132.obsm['X_umap_transfer'] = X132_umap
print('  UMAP done.')

# ══════════════════════════════════════════════════════════════════════════
#  Step 5: ABMIL 推断 GSE132465（冻结权重）
# ══════════════════════════════════════════════════════════════════════════
print('\nStep 5: ABMIL inference on GSE132465 (frozen weights) ...')
pids132     = adata132.obs['Patient'].values
unique132   = adata132.obs['Patient'].unique()
bags132     = {pid: {'idx': np.where(pids132 == pid)[0]} for pid in unique132}
per_patient = {pid: len(v['idx']) for pid, v in bags132.items()}
print(f'  {len(unique132)} patients, median {np.median(list(per_patient.values())):.0f} cells/patient')

final_model.eval()
att132 = np.zeros(adata132.n_obs, dtype=np.float32)
pred_probs = {}

with torch.no_grad():
    for pid in unique132:
        idx = bags132[pid]['idx']
        x   = torch.tensor(X132_pca[idx], dtype=torch.float32)
        logit, A = final_model(x, return_attention=True)
        att132[idx] = A.numpy()
        pred_probs[pid] = torch.sigmoid(logit).item()

adata132.obs['attention_score'] = att132
print(f'  Attention: min={att132.min():.6f}  max={att132.max():.6f}  mean={att132.mean():.6f}')

# 模型预测（MMRd 概率，仅供参考）
pred_df = pd.DataFrame({'patient': list(pred_probs.keys()),
                         'prob_MMRd': list(pred_probs.values()),
                         'n_cells': [per_patient[p] for p in pred_probs.keys()]})
pred_df = pred_df.sort_values('prob_MMRd', ascending=False)
print('\n  GSE132465 predicted MMRd probability (top 10):')
print(pred_df.head(10).to_string(index=False))
pred_df.to_csv(f'{RES_DIR}/gse132465_predictions.csv', index=False)

# ── 按亚型统计 attention ──────────────────────────────────────────────────
adata132.obs['subtype_mapped'] = (adata132.obs['unified_cell_subtype']
                                    .map(GSE132_SUBTYPE_MAP)
                                    .fillna('Other'))
att132_by_ct = adata132.obs.groupby('subtype_mapped')['attention_score'].agg(
    ['mean', 'median', 'std', 'count']).reset_index()
att132_by_ct.columns = ['subtype', 'mean', 'median', 'std', 'n']
att132_by_ct = att132_by_ct[att132_by_ct['n'] >= 50].sort_values('mean', ascending=False)
print(f'\n  GSE132465 attention by subtype:\n{att132_by_ct.to_string(index=False)}')
att132_by_ct.to_csv(f'{RES_DIR}/gse132465_att_by_subtype.csv', index=False)

# ══════════════════════════════════════════════════════════════════════════
#  Panel C: Attention UMAP（2 subplots: attention score + subtype）
# ══════════════════════════════════════════════════════════════════════════
print('\nPlotting Panel C ...')

fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
umap_xy = X132_umap

# ─── 左: Attention Score ─────────────────────────────────────────────────
ax = axes[0]
# 先画低分（背景），再画高分（前景）
sort_order = np.argsort(att132)
sc1 = ax.scatter(umap_xy[sort_order, 0], umap_xy[sort_order, 1],
                  c=att132[sort_order],
                  cmap='plasma', s=1.5, linewidths=0, alpha=0.7,
                  vmin=np.percentile(att132, 2),
                  vmax=np.percentile(att132, 99))
cb = plt.colorbar(sc1, ax=ax, shrink=0.75, aspect=20, pad=0.02)
cb.set_label('Attention Score', fontsize=9)
cb.ax.tick_params(labelsize=7.5)
ax.set_title('Attention Score\nGSE132465 · Immune Cells (Tumor)',
             fontsize=11, fontweight='bold')
ax.set_xlabel('UMAP 1', fontsize=9)
ax.set_ylabel('UMAP 2', fontsize=9)
ax.tick_params(labelsize=8)
ax.set_aspect('equal')

# ─── 右: Cell Subtype ────────────────────────────────────────────────────
ax = axes[1]
subtypes = adata132.obs['subtype_mapped'].values
uniq_sub = [s for s in ['Macro', 'TCD8', 'TCD4', 'Treg', 'NK', 'DC', 'Tgd', 'Other']
             if s in set(subtypes)]

for sub in uniq_sub:
    mask = subtypes == sub
    if mask.sum() == 0:
        continue
    color = SUBTYPE_COLORS.get(sub, '#DDDDDD')
    size  = 4.0 if sub in PRE_DRIVER else 2.0
    alpha = 0.85 if sub in PRE_DRIVER else 0.55
    zord  = 5 if sub in PRE_DRIVER else 2
    label = f'{sub} (Pre-driver)' if sub in PRE_DRIVER else sub
    ax.scatter(umap_xy[mask, 0], umap_xy[mask, 1],
               c=color, s=size, linewidths=0, alpha=alpha,
               zorder=zord, label=label)

ax.set_title('Cell Subtype\nGSE132465 · Immune Cells (Tumor)',
             fontsize=11, fontweight='bold')
ax.set_xlabel('UMAP 1', fontsize=9)
ax.set_ylabel('UMAP 2', fontsize=9)
ax.tick_params(labelsize=8)
ax.set_aspect('equal')
ax.legend(fontsize=7.5, markerscale=3, loc='upper left',
          framealpha=0.8, edgecolor='lightgray',
          handlelength=1.0, handletextpad=0.5)

fig.suptitle('Weight Transfer Validation — GSE132465 (Independent CRC Cohort)',
             fontsize=12, fontweight='bold', y=1.01)

plt.tight_layout()
fig.savefig(f'{FIG_DIR}/figure3C_attention_umap.png', dpi=200, bbox_inches='tight')
fig.savefig(f'{FIG_DIR}/figure3C_attention_umap.pdf', bbox_inches='tight')
plt.close(fig)
print('  Saved figure3C_attention_umap.png/.pdf')

# ══════════════════════════════════════════════════════════════════════════
#  Panel D: Subtype Attention Ranking (Discovery vs Validation)
# ══════════════════════════════════════════════════════════════════════════
print('\nPlotting Panel D ...')

# 对齐两数据集的亚型名称
sub178 = att178_by_ct.set_index('subtype')[['mean', 'std']].rename(
    columns={'mean': 'mean_178', 'std': 'std_178'})
sub132 = att132_by_ct.set_index('subtype')[['mean', 'std']].rename(
    columns={'mean': 'mean_132', 'std': 'std_132'})
comp = sub178.join(sub132, how='outer').fillna(0)

# 按 GSE178341（发现集）排序
comp = comp.sort_values('mean_178', ascending=True).reset_index()

n_sub = len(comp)
fig, ax = plt.subplots(figsize=(7, max(4.5, n_sub * 0.55 + 1.2)))

y = np.arange(n_sub)
offset = 0.16  # 两数据集点的垂直偏移

for i, row in comp.iterrows():
    ct = row['subtype']
    is_pd = ct in PRE_DRIVER
    c178 = SUBTYPE_COLORS.get(ct, '#888888') if is_pd else '#3C5488'
    c132 = SUBTYPE_COLORS.get(ct, '#888888') if is_pd else '#91D1C2'

    # GSE178341（实心圆）
    if row['mean_178'] > 0:
        ax.plot([0, row['mean_178']], [i + offset, i + offset],
                color=c178, lw=1.5, alpha=0.8, zorder=2)
        ax.scatter(row['mean_178'], i + offset,
                   s=80 if is_pd else 50, c=c178,
                   zorder=3, edgecolors='white', linewidths=0.5)

    # GSE132465（空心菱形）
    if row['mean_132'] > 0:
        ax.plot([0, row['mean_132']], [i - offset, i - offset],
                color=c132, lw=1.5, alpha=0.8, ls='--', zorder=2)
        ax.scatter(row['mean_132'], i - offset,
                   s=80 if is_pd else 50, c=c132, marker='D',
                   zorder=3, edgecolors='white', linewidths=0.5)

ax.axvline(0, color='black', lw=0.6)
ax.set_yticks(y)
ax.set_yticklabels(comp['subtype'], fontsize=10)
for tick, ct in zip(ax.get_yticklabels(), comp['subtype']):
    if ct in PRE_DRIVER:
        tick.set_color(SUBTYPE_COLORS.get(ct, '#E64B35'))
        tick.set_fontweight('bold')

# Pre-driver 右侧标注
for i, row in comp.iterrows():
    if row['subtype'] in PRE_DRIVER:
        ax.text(1.02, i, 'Pre-driver', transform=ax.get_yaxis_transform(),
                va='center', ha='left', fontsize=7.5,
                color=SUBTYPE_COLORS.get(row['subtype'], '#E64B35'),
                fontstyle='italic')

ax.set_xlabel('Mean Attention Score per Cell', fontsize=11)
ax.set_title('Immune Subtype Attention Weight Ranking\nDiscovery (GSE178341) vs Validation (GSE132465)',
             fontsize=12, fontweight='bold', pad=12)
ax.set_clip_on(False)

# 图例
legend_patches = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#3C5488',
               markersize=8, label='GSE178341 (Discovery · MMRp)'),
    plt.Line2D([0], [0], marker='D', color='w', markerfacecolor='#91D1C2',
               markersize=7, label='GSE132465 (Validation · CRC Tumor)'),
    mpatches.Patch(color='#E64B35', label='Pre-driver: Macro'),
    mpatches.Patch(color='#3C5488', label='Pre-driver: TCD8'),
]
ax.legend(handles=legend_patches, fontsize=8.5, loc='lower right',
          framealpha=0.85, edgecolor='lightgray')

plt.tight_layout()
fig.savefig(f'{FIG_DIR}/figure3D_subtype_weight_ranking.png', dpi=200, bbox_inches='tight')
fig.savefig(f'{FIG_DIR}/figure3D_subtype_weight_ranking.pdf', bbox_inches='tight')
plt.close(fig)
print('  Saved figure3D_subtype_weight_ranking.png/.pdf')

# ══════════════════════════════════════════════════════════════════════════
#  Summary
# ══════════════════════════════════════════════════════════════════════════
print('\n' + '=' * 65)
print('SUMMARY')
print(f'\nCV AUC on GSE178341 (common HVG PCA): {np.mean(fold_aucs):.4f} ± {np.std(fold_aucs):.4f}')
print(f'\nGSE132465 attention by subtype (sorted):')
print(att132_by_ct[['subtype', 'mean', 'n']].to_string(index=False))
print('\nPre-driver subtypes in GSE132465:')
for _, r in att132_by_ct[att132_by_ct['subtype'].isin(PRE_DRIVER)].iterrows():
    print(f'  {r["subtype"]:8s}: mean={r["mean"]:.6f}  n={r["n"]:,}')
print('\nDone.')
