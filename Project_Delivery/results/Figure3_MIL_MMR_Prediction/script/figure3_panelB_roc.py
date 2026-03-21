#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 3 Panel B — ABMIL Performance: ROC + PR Curves
GSE178341 immune-only bags · 5-fold CV out-of-fold predictions
"""

import scanpy as sc
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.optim import Adam
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (roc_auc_score, roc_curve,
                              precision_recall_curve, average_precision_score)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
os.makedirs(FIG_DIR, exist_ok=True)

# ── 超参数（与 phase3_mil_immune.py 完全一致）─────────────────────────────
N_PCA      = 50
MAX_CELLS  = 1500
HIDDEN_DIM = 256
ATT_DIM    = 128
DROPOUT    = 0.25
LR         = 1e-4
WD         = 1e-5
EPOCHS     = 100
N_FOLDS    = 5
PATIENCE   = 20

# ══════════════════════════════════════════════════════════════════════════
#  Step 1: 加载数据，重建 bags
# ══════════════════════════════════════════════════════════════════════════
print('Loading data ...')
adata = sc.read_h5ad('processed_data/GSE178341_immune_mil.h5ad')
print(f'  {adata.n_obs:,} cells, {adata.obs["PID"].nunique()} patients')

pca_feat = adata.obsm['X_pca'][:, :N_PCA].astype('float32')
pids     = adata.obs['PID'].values
mmr_lbl  = (adata.obs['MMRStatus'] == 'MMRd').astype(int).values

unique_pids = adata.obs['PID'].unique()
patient_bags = {}
for pid in unique_pids:
    idx = np.where(pids == pid)[0]
    patient_bags[pid] = {'idx': idx, 'label': int(mmr_lbl[idx[0]])}

all_pids   = np.array(list(patient_bags.keys()))
all_labels = np.array([patient_bags[p]['label'] for p in all_pids])
print(f'  Bags: {len(all_pids)}  (MMRd={all_labels.sum()}, MMRp={(1-all_labels).sum()})')

# ── Model ─────────────────────────────────────────────────────────────────
class GatedABMIL(nn.Module):
    def __init__(self):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(N_PCA, HIDDEN_DIM), nn.ReLU(), nn.Dropout(DROPOUT),
            nn.Linear(HIDDEN_DIM, HIDDEN_DIM), nn.ReLU(),
        )
        self.att_V = nn.Sequential(nn.Linear(HIDDEN_DIM, ATT_DIM), nn.Tanh())
        self.att_U = nn.Sequential(nn.Linear(HIDDEN_DIM, ATT_DIM), nn.Sigmoid())
        self.att_w = nn.Linear(ATT_DIM, 1, bias=False)
        self.classifier = nn.Sequential(nn.Dropout(DROPOUT), nn.Linear(HIDDEN_DIM, 1))

    def forward(self, x):
        H = self.encoder(x)
        A = F.softmax(self.att_w(self.att_V(H) * self.att_U(H)), dim=0)
        return self.classifier((A * H).sum(dim=0, keepdim=True)).squeeze()


def train_epoch(model, opt, pids_list):
    model.train()
    order = pids_list.copy()
    np.random.shuffle(order)
    for pid in order:
        idx = patient_bags[pid]['idx']
        if len(idx) > MAX_CELLS:
            idx = np.random.choice(idx, MAX_CELLS, replace=False)
        x   = torch.tensor(pca_feat[idx])
        lbl = torch.tensor(patient_bags[pid]['label'], dtype=torch.float32)
        opt.zero_grad()
        F.binary_cross_entropy_with_logits(model(x), lbl).backward()
        opt.step()


@torch.no_grad()
def predict(model, pids_list):
    model.eval()
    probs, labels = [], []
    for pid in pids_list:
        idx = patient_bags[pid]['idx']
        probs.append(torch.sigmoid(model(torch.tensor(pca_feat[idx]))).item())
        labels.append(patient_bags[pid]['label'])
    return np.array(probs), np.array(labels)


# ══════════════════════════════════════════════════════════════════════════
#  Step 2: 5-Fold CV — 收集 out-of-fold 预测
# ══════════════════════════════════════════════════════════════════════════
print('\n5-Fold CV (collecting out-of-fold predictions) ...')
skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=42)

oof_probs   = np.zeros(len(all_pids))
oof_labels  = np.zeros(len(all_pids), dtype=int)
fold_roc    = []   # (fpr, tpr, auc) per fold
fold_pr     = []   # (recall, precision, ap) per fold

for fold, (tr_idx, va_idx) in enumerate(skf.split(all_pids, all_labels)):
    pids_tr = all_pids[tr_idx].tolist()
    pids_va = all_pids[va_idx].tolist()

    model = GatedABMIL()
    opt   = Adam(model.parameters(), lr=LR, weight_decay=WD)
    best_auc, best_state, patience_cnt = 0.0, None, 0

    for epoch in range(EPOCHS):
        train_epoch(model, opt, pids_tr)
        probs_va, lbl_va = predict(model, pids_va)
        if len(np.unique(lbl_va)) > 1:
            auc = roc_auc_score(lbl_va, probs_va)
        else:
            auc = 0.5
        if auc > best_auc:
            best_auc   = auc
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
            patience_cnt = 0
        else:
            patience_cnt += 1
        if patience_cnt >= PATIENCE:
            break

    model.load_state_dict(best_state)
    probs_va, lbl_va = predict(model, pids_va)

    # store OOF
    oof_probs[va_idx]  = probs_va
    oof_labels[va_idx] = lbl_va

    # fold-level curves
    if len(np.unique(lbl_va)) > 1:
        fpr, tpr, _ = roc_curve(lbl_va, probs_va)
        fold_roc.append((fpr, tpr, roc_auc_score(lbl_va, probs_va)))
        prec, rec, _ = precision_recall_curve(lbl_va, probs_va)
        fold_pr.append((rec, prec, average_precision_score(lbl_va, probs_va)))

    print(f'  Fold {fold+1}: AUC={best_auc:.4f}  n_val={len(pids_va)}')

# 全局 OOF 指标
oof_auc = roc_auc_score(oof_labels, oof_probs)
oof_ap  = average_precision_score(oof_labels, oof_probs)
print(f'\nOOF AUC = {oof_auc:.4f}   OOF AP (AUPRC) = {oof_ap:.4f}')

# ══════════════════════════════════════════════════════════════════════════
#  Step 3: Plot Panel B — ROC + PR
# ══════════════════════════════════════════════════════════════════════════
print('\nPlotting Panel B ...')

COLOR_MAIN  = '#3C5488'   # 主线（深蓝）
COLOR_FOLD  = '#8491B4'   # 折叠线（浅蓝）
COLOR_DIAG  = '#AAAAAA'   # 对角线

fig, axes = plt.subplots(1, 2, figsize=(8.5, 4.0))

# ─── 左图：ROC ───────────────────────────────────────────────────────────
ax = axes[0]

# 各折细线
for fpr, tpr, auc_f in fold_roc:
    ax.plot(fpr, tpr, color=COLOR_FOLD, lw=0.9, alpha=0.45)

# 全局 OOF ROC
fpr_oof, tpr_oof, _ = roc_curve(oof_labels, oof_probs)
ax.plot(fpr_oof, tpr_oof, color=COLOR_MAIN, lw=2.2,
        label=f'OOF  AUC = {oof_auc:.3f}')

# 对角线
ax.plot([0, 1], [0, 1], ls='--', color=COLOR_DIAG, lw=1.0, label='Random')

# 阴影：fold range（用 interpolation）
from scipy.interpolate import interp1d
common_fpr = np.linspace(0, 1, 300)
tpr_mat = []
for fpr_f, tpr_f, _ in fold_roc:
    f_interp = interp1d(fpr_f, tpr_f, kind='linear',
                        bounds_error=False, fill_value=(0, 1))
    tpr_mat.append(f_interp(common_fpr))
tpr_mat = np.array(tpr_mat)
ax.fill_between(common_fpr,
                tpr_mat.min(axis=0), tpr_mat.max(axis=0),
                alpha=0.12, color=COLOR_MAIN, label='Fold range')

ax.set_xlim(-0.02, 1.02)
ax.set_ylim(-0.02, 1.05)
ax.set_xlabel('False Positive Rate', fontsize=11)
ax.set_ylabel('True Positive Rate', fontsize=11)
ax.set_title('ROC Curve', fontsize=12, fontweight='bold')
ax.legend(fontsize=9, loc='lower right', framealpha=0.85, edgecolor='lightgray')
ax.set_aspect('equal')

# ─── 右图：PR ────────────────────────────────────────────────────────────
ax = axes[1]

# 各折细线
for rec_f, prec_f, ap_f in fold_pr:
    ax.plot(rec_f, prec_f, color=COLOR_FOLD, lw=0.9, alpha=0.45)

# 全局 OOF PR
prec_oof, rec_oof, _ = precision_recall_curve(oof_labels, oof_probs)
ax.plot(rec_oof, prec_oof, color=COLOR_MAIN, lw=2.2,
        label=f'OOF  AP = {oof_ap:.3f}')

# 基线（MMRd 比例）
baseline = oof_labels.mean()
ax.axhline(baseline, ls='--', color=COLOR_DIAG, lw=1.0,
           label=f'Baseline (prevalence = {baseline:.2f})')

# 阴影：fold range
common_rec = np.linspace(0, 1, 300)
prec_mat = []
for rec_f, prec_f, _ in fold_pr:
    # PR 曲线从右到左，需要反转排序
    sort_idx = np.argsort(rec_f)
    f_interp = interp1d(rec_f[sort_idx], prec_f[sort_idx], kind='linear',
                        bounds_error=False, fill_value=(prec_f[-1], prec_f[0]))
    prec_mat.append(f_interp(common_rec))
prec_mat = np.array(prec_mat)
ax.fill_between(common_rec,
                prec_mat.min(axis=0), prec_mat.max(axis=0),
                alpha=0.12, color=COLOR_MAIN, label='Fold range')

ax.set_xlim(-0.02, 1.02)
ax.set_ylim(0.0, 1.05)
ax.set_xlabel('Recall', fontsize=11)
ax.set_ylabel('Precision', fontsize=11)
ax.set_title('Precision–Recall Curve', fontsize=12, fontweight='bold')
ax.legend(fontsize=9, loc='upper right', framealpha=0.85, edgecolor='lightgray')
ax.set_aspect('equal')

# ─── 总标题 ──────────────────────────────────────────────────────────────
fig.suptitle('ABMIL Performance · GSE178341 Immune Microenvironment\n'
             'MMRp vs MMRd · 5-Fold Cross-Validation',
             fontsize=11, fontweight='bold', y=1.03)

plt.tight_layout()
fig.savefig(f'{FIG_DIR}/figure3B_panelB_roc_pr.png', dpi=200, bbox_inches='tight')
fig.savefig(f'{FIG_DIR}/figure3B_panelB_roc_pr.pdf', bbox_inches='tight')
plt.close(fig)
print(f'Saved: figures/phase3/figure3B_panelB_roc_pr.png/.pdf')
print('Done.')
