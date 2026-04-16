#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
External validation using GSE146771 (CRC leukocytes, 43,817 cells)
Strategy:
  1. Convert GSE146771 log1p(TPM) -> approximate counts
  2. Subset to common genes (8,780)
  3. Embed via scVI encoder (encoder is batch-agnostic)
  4. KNN classify against annotated training set
  5. Compare predictions to GSE146771 ground truth
"""

import scanpy as sc
import scvi
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import (classification_report, confusion_matrix,
                             accuracy_score, f1_score)
import os, sys, io
from datetime import datetime

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

OUT_DIR = 'validation_results'
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 70)
print("External Validation: GSE146771")
print("=" * 70)

# ─────────────────────────────────────────────────────────────────────────────
# 1. Load data
# ─────────────────────────────────────────────────────────────────────────────
print("\n[1/7] Loading datasets ...")

ref   = sc.read_h5ad('processed_data/merged_annotated.h5ad')
query = sc.read_h5ad('processed_data/GSE146771_tpm.h5ad')

print(f"  Reference (training): {ref.n_obs:,} cells x {ref.n_vars:,} genes")
print(f"  Query (GSE146771):    {query.n_obs:,} cells x {query.n_vars:,} genes")

# ─────────────────────────────────────────────────────────────────────────────
# 2. Prepare GSE146771: log1p(TPM) -> approximate counts
# ─────────────────────────────────────────────────────────────────────────────
print("\n[2/7] Converting GSE146771 log1p(TPM) -> approximate counts ...")

# Convert: count_approx = expm1(X) * raw_nUMI / 1e6
raw_numi = query.obs['raw.nUMI'].values
X_tpm    = query.X.toarray() if hasattr(query.X, 'toarray') else query.X.copy()
X_tpm    = X_tpm.astype(np.float32)

# Back-transform: TPM = expm1(log1p(TPM))
X_tpm_raw = np.expm1(X_tpm)

# Approx counts per cell
X_counts = X_tpm_raw * raw_numi[:, None] / 1e6
X_counts = np.round(X_counts).astype(np.float32)
X_counts = np.clip(X_counts, 0, None)

import scipy.sparse as sp
query_counts = query.copy()
query_counts.X = sp.csr_matrix(X_counts)

# QC check
approx_umis = np.array(query_counts.X.sum(axis=1)).flatten()
print(f"  Approx UMI per cell: median={np.median(approx_umis):.0f}  "
      f"mean={approx_umis.mean():.0f}  max={approx_umis.max():.0f}")

# ─────────────────────────────────────────────────────────────────────────────
# 3. Align genes: subset to common genes in same order as training
# ─────────────────────────────────────────────────────────────────────────────
print("\n[3/7] Aligning gene space ...")

common_genes = ref.var_names.intersection(query_counts.var_names)
print(f"  Common genes: {len(common_genes):,} / {ref.n_vars:,} training genes "
      f"({len(common_genes)/ref.n_vars*100:.1f}%)")

# Build query with ALL 9,009 training genes (pad missing 229 with zeros)
import scipy.sparse as sp_build
import anndata as ad

missing_genes = ref.var_names.difference(query_counts.var_names)
print(f"  Missing genes (set to 0): {len(missing_genes):,}")

# Start from common-gene subset, then add zero columns for missing genes
query_common = query_counts[:, common_genes].copy()

# Build zero block for missing genes
n_cells    = query_common.n_obs
n_missing  = len(missing_genes)
zero_block = sp_build.csr_matrix((n_cells, n_missing), dtype=np.float32)

# Create AnnData for missing genes
zero_adata = ad.AnnData(
    X   = zero_block,
    obs = query_common.obs,
    var = pd.DataFrame(index=missing_genes),
)

# Concatenate: common genes + zero missing genes, then reorder to match training
query_full = ad.concat([query_common, zero_adata], axis=1)
query_full = query_full[:, ref.var_names].copy()   # reorder to exact training order

print(f"  Query aligned to {query_full.n_vars:,} genes (training order)")

query_sub = query_full
ref_sub   = ref   # use full reference (9,009 genes)

# Add required obs columns for scVI setup
query_sub.obs['batch']   = 'GSE132257'   # CRC reference batch (encoder is batch-agnostic)
query_sub.obs['disease'] = 'CRC'

# ─────────────────────────────────────────────────────────────────────────────
# 4. Embed GSE146771 via scVI encoder
# ─────────────────────────────────────────────────────────────────────────────
print("\n[4/7] Embedding GSE146771 via scVI encoder ...")

# Re-setup with common genes only (need consistent var_names)
scvi.model.SCVI.setup_anndata(
    query_sub,
    batch_key='batch',
    categorical_covariate_keys=['disease']
)

model = scvi.model.SCVI.load('scvi_model_optimized', adata=ref_sub)
print(f"  Model loaded.")

query_latent = model.get_latent_representation(adata=query_sub)
print(f"  Query embeddings: {query_latent.shape}")

# ─────────────────────────────────────────────────────────────────────────────
# 5. KNN classification
# ─────────────────────────────────────────────────────────────────────────────
print("\n[5/7] KNN classification ...")

# Reference embeddings + labels (immune cells only for fair comparison)
immune_types = {'T cells', 'B cells', 'NK cells', 'Myeloids'}
ref_immune_mask = ref.obs['unified_cell_type'].isin(immune_types)
X_ref   = ref.obsm['X_scvi'][ref_immune_mask.values]
y_ref   = ref.obs.loc[ref_immune_mask, 'unified_cell_type'].values
print(f"  Reference immune cells: {len(X_ref):,}  (types: {np.unique(y_ref).tolist()})")

knn = KNeighborsClassifier(n_neighbors=15, metric='euclidean', n_jobs=-1)
knn.fit(X_ref, y_ref)

y_pred = knn.predict(query_latent)
y_prob = knn.predict_proba(query_latent)
confidence = y_prob.max(axis=1)

# Ground-truth mapping: GSE146771 labels -> our label space
gt_map = {
    'B cell':       'B cells',
    'CD4 T cell':   'T cells',
    'CD8 T cell':   'T cells',
    'ILC':          'NK cells',
    'Myeloid cell': 'Myeloids',
}
y_true = query.obs['Global_Cluster'].map(gt_map).values
print(f"  GT mapping applied: {gt_map}")

# ─────────────────────────────────────────────────────────────────────────────
# 6. Metrics
# ─────────────────────────────────────────────────────────────────────────────
print("\n[6/7] Computing metrics ...")

label_order = ['T cells', 'B cells', 'NK cells', 'Myeloids']

acc    = accuracy_score(y_true, y_pred)
f1_mac = f1_score(y_true, y_pred, average='macro',    labels=label_order, zero_division=0)
f1_wt  = f1_score(y_true, y_pred, average='weighted', labels=label_order, zero_division=0)

print(f"\n  Overall accuracy : {acc:.4f} ({acc*100:.1f}%)")
print(f"  F1 macro         : {f1_mac:.4f}")
print(f"  F1 weighted      : {f1_wt:.4f}")
print(f"  Mean confidence  : {confidence.mean():.4f}")
print()
print("  Per-class report:")
print(classification_report(y_true, y_pred, labels=label_order, zero_division=0))

# Per-original-label breakdown
print("  Accuracy by GSE146771 original label:")
for orig_label, mapped in gt_map.items():
    mask = (query.obs['Global_Cluster'] == orig_label).values
    if mask.sum() == 0:
        continue
    sub_acc = accuracy_score(y_true[mask], y_pred[mask])
    print(f"    {orig_label:20s} -> {mapped:12s}: "
          f"n={mask.sum():,}  acc={sub_acc:.3f}")

# Sub-cluster accuracy (fine-grained)
print("\n  Top-10 Sub_Cluster accuracy:")
sub_df = pd.DataFrame({
    'Sub_Cluster': query.obs['Sub_Cluster'].values,
    'y_true': y_true,
    'y_pred': y_pred,
})
sub_acc_df = (sub_df.groupby('Sub_Cluster')
              .apply(lambda g: pd.Series({
                  'n': len(g),
                  'acc': accuracy_score(g['y_true'], g['y_pred'])
              }))
              .sort_values('n', ascending=False)
              .head(15))
print(sub_acc_df.to_string())

# ─────────────────────────────────────────────────────────────────────────────
# 7. Plots
# ─────────────────────────────────────────────────────────────────────────────
print("\n[7/7] Generating plots ...")

fig = plt.figure(figsize=(20, 16))
gs  = fig.add_gridspec(2, 3, hspace=0.4, wspace=0.35)

# ── Panel A: Confusion matrix (normalised) ────────────────────────────────
ax = fig.add_subplot(gs[0, 0])
cm = confusion_matrix(y_true, y_pred, labels=label_order, normalize='true')
sns.heatmap(cm, annot=True, fmt='.2f', cmap='Blues',
            xticklabels=label_order, yticklabels=label_order, ax=ax,
            cbar_kws={'shrink': 0.8}, linewidths=0.5)
ax.set_title(f'Confusion Matrix\n(normalised by row)\nOverall acc = {acc:.3f}',
             fontsize=11, fontweight='bold')
ax.set_xlabel('Predicted', fontsize=10)
ax.set_ylabel('True (GSE146771)', fontsize=10)
ax.tick_params(axis='x', rotation=30, labelsize=8)
ax.tick_params(axis='y', rotation=0,  labelsize=8)

# ── Panel B: Per-original-label accuracy ────────────────────────────────────
ax = fig.add_subplot(gs[0, 1])
orig_labels = list(gt_map.keys())
orig_accs   = []
orig_ns     = []
for lbl in orig_labels:
    mask = (query.obs['Global_Cluster'] == lbl).values
    orig_accs.append(accuracy_score(y_true[mask], y_pred[mask]))
    orig_ns.append(mask.sum())
colors_bar = ['#4DBBD5' if a >= 0.8 else '#E64B35' if a < 0.6 else '#F39B7F'
               for a in orig_accs]
bars = ax.barh(orig_labels, orig_accs, color=colors_bar, edgecolor='white')
ax.axvline(0.8, color='#E64B35', linestyle='--', linewidth=1.2, label='80% threshold')
ax.set_xlim(0, 1.05)
ax.set_xlabel('Accuracy', fontsize=10)
ax.set_title('Accuracy by Original\nGSE146771 Cell Type', fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
for bar, acc_val, n in zip(bars, orig_accs, orig_ns):
    ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
            f'{acc_val:.2f}  (n={n:,})', va='center', fontsize=9)

# ── Panel C: Confidence distribution ────────────────────────────────────────
ax = fig.add_subplot(gs[0, 2])
correct = (y_true == y_pred)
ax.hist(confidence[correct],  bins=40, alpha=0.7, color='#4DBBD5',
        label=f'Correct ({correct.sum():,})')
ax.hist(confidence[~correct], bins=40, alpha=0.7, color='#E64B35',
        label=f'Wrong   ({(~correct).sum():,})')
ax.set_xlabel('KNN Confidence', fontsize=10)
ax.set_ylabel('Cells', fontsize=10)
ax.set_title('KNN Prediction Confidence\n(correct vs wrong)',
             fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.axvline(0.5, color='grey', linestyle='--', linewidth=1)

# ── Panel D: Sub-cluster accuracy heatmap ────────────────────────────────────
ax = fig.add_subplot(gs[1, :2])
sub_full = (sub_df.groupby(['Sub_Cluster', 'y_pred'])
            .size()
            .unstack(fill_value=0))
# Normalise by row
sub_full_norm = sub_full.div(sub_full.sum(axis=1), axis=0)
# Keep top 20 by total cells
top_subs = sub_df['Sub_Cluster'].value_counts().head(20).index
sub_plot = sub_full_norm.loc[sub_full_norm.index.isin(top_subs)]
sub_plot = sub_plot.reindex(columns=label_order, fill_value=0)
# Sort by dominant type
sub_plot = sub_plot.sort_values(label_order[0], ascending=False)
sns.heatmap(sub_plot, annot=True, fmt='.2f', cmap='YlOrRd',
            ax=ax, linewidths=0.3, cbar_kws={'shrink': 0.6},
            vmin=0, vmax=1)
ax.set_title('Sub-cluster Prediction Distribution\n(top 20 sub-clusters, normalised)',
             fontsize=11, fontweight='bold')
ax.set_xlabel('Our Predicted Label', fontsize=10)
ax.set_ylabel('GSE146771 Sub-Cluster', fontsize=10)
ax.tick_params(axis='y', labelsize=8)
ax.tick_params(axis='x', rotation=30, labelsize=9)

# ── Panel E: Summary text ────────────────────────────────────────────────────
ax = fig.add_subplot(gs[1, 2])
ax.axis('off')
summary_text = (
    f"External Validation Summary\n"
    f"{'─'*32}\n"
    f"Query:  GSE146771 (CRC leukocytes)\n"
    f"        {query.n_obs:,} cells, 5 types\n\n"
    f"Reference: 5-batch training set\n"
    f"           {ref_immune_mask.sum():,} immune cells\n\n"
    f"Method: scVI encoder + KNN (k=15)\n"
    f"Genes:  {len(common_genes):,} / {ref.n_vars:,} ({len(common_genes)/ref.n_vars*100:.1f}%)\n\n"
    f"{'─'*32}\n"
    f"Overall accuracy:   {acc*100:.1f}%\n"
    f"F1 macro:           {f1_mac:.3f}\n"
    f"F1 weighted:        {f1_wt:.3f}\n"
    f"Mean confidence:    {confidence.mean():.3f}\n\n"
    f"Per-type accuracy:\n"
)
for lbl, acc_val, n in zip(orig_labels, orig_accs, orig_ns):
    summary_text += f"  {lbl:<20s}: {acc_val*100:.1f}%\n"

ax.text(0.05, 0.98, summary_text, transform=ax.transAxes,
        fontsize=10, verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='#f0f4f8', alpha=0.8))

fig.suptitle('External Validation: GSE146771 (CRC Leukocytes)',
             fontsize=15, fontweight='bold', y=1.01)

out_path = f'{OUT_DIR}/external_validation.png'
plt.savefig(out_path, dpi=200, bbox_inches='tight')
plt.close()
print(f"  Saved: {out_path}")

# ── Save result table ─────────────────────────────────────────────────────────
result_df = pd.DataFrame({
    'cell_id':         query.obs.index,
    'gse_true_label':  query.obs['Global_Cluster'].values,
    'gse_sub_cluster': query.obs['Sub_Cluster'].values,
    'mapped_true':     y_true,
    'our_predicted':   y_pred,
    'confidence':      confidence,
    'correct':         (y_true == y_pred),
})
result_df.to_csv(f'{OUT_DIR}/validation_predictions.csv', index=False)
print(f"  Saved: {OUT_DIR}/validation_predictions.csv")

print("\n" + "=" * 70)
print(f"Validation COMPLETE — {datetime.now().strftime('%H:%M:%S')}")
print(f"  Overall accuracy : {acc*100:.1f}%")
print(f"  F1 macro         : {f1_mac:.3f}")
print("=" * 70)
