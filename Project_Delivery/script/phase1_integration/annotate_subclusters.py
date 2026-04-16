#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reverse KNN: use GSE146771 Sub_Cluster (38 types) as reference,
annotate immune cells in the 5-batch training set.
"""

import scanpy as sc
import scvi
import pandas as pd
import numpy as np
import anndata as ad
import scipy.sparse as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import classification_report, accuracy_score
import os, sys, io
from datetime import datetime

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

OUT_DIR = 'annotation_results'
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 70)
print("Sub-cluster Annotation via GSE146771 (reverse KNN)")
print("=" * 70)

# ─────────────────────────────────────────────────────────────────────────────
# 1. Load data
# ─────────────────────────────────────────────────────────────────────────────
print("\n[1/6] Loading data ...")
ref   = sc.read_h5ad('processed_data/merged_annotated.h5ad')   # 5-batch, has X_scvi
query = sc.read_h5ad('processed_data/GSE146771_tpm.h5ad')      # reference for sub-clusters

print(f"  Training set: {ref.n_obs:,} cells, X_scvi={ref.obsm['X_scvi'].shape}")
print(f"  GSE146771:    {query.n_obs:,} cells, {query.obs['Sub_Cluster'].nunique()} sub-clusters")
print(f"\n  All 38 sub-clusters:")
for sc_name, n in query.obs['Sub_Cluster'].value_counts().items():
    print(f"    {sc_name:35s}: {n:,}")

# ─────────────────────────────────────────────────────────────────────────────
# 2. Re-embed GSE146771 via scVI encoder
# ─────────────────────────────────────────────────────────────────────────────
print("\n[2/6] Re-embedding GSE146771 via scVI encoder ...")

# log1p(TPM) -> approximate counts
raw_numi = query.obs['raw.nUMI'].values
X_tpm    = query.X.toarray().astype(np.float32)
X_counts = np.round(np.clip(np.expm1(X_tpm) * raw_numi[:, None] / 1e6, 0, None)).astype(np.float32)

# Align to training gene space (pad missing genes with 0)
common_genes = ref.var_names.intersection(query.var_names)
missing_genes = ref.var_names.difference(query.var_names)
print(f"  Common genes: {len(common_genes):,}  |  Padded with 0: {len(missing_genes):,}")

query_common = ad.AnnData(
    X   = sp.csr_matrix(X_counts[:, query.var_names.get_indexer(common_genes)]),
    obs = query.obs.copy(),
    var = pd.DataFrame(index=common_genes),
)
zero_block  = sp.csr_matrix((query.n_obs, len(missing_genes)), dtype=np.float32)
zero_adata  = ad.AnnData(
    X   = zero_block,
    obs = query.obs.copy(),
    var = pd.DataFrame(index=missing_genes),
)
query_full = ad.concat([query_common, zero_adata], axis=1)[:, ref.var_names].copy()
query_full.obs['batch']   = 'GSE132257'
query_full.obs['disease'] = 'CRC'

scvi.model.SCVI.setup_anndata(query_full, batch_key='batch',
                               categorical_covariate_keys=['disease'])
model = scvi.model.SCVI.load('scvi_model_optimized', adata=ref)
gse_latent = model.get_latent_representation(adata=query_full)
print(f"  GSE146771 embedding: {gse_latent.shape}")

# ─────────────────────────────────────────────────────────────────────────────
# 3. Train KNN on GSE146771 sub-clusters (reference)
# ─────────────────────────────────────────────────────────────────────────────
print("\n[3/6] Training KNN on GSE146771 sub-clusters ...")

X_gse   = gse_latent                              # (43817, 30)  reference
y_gse   = query.obs['Sub_Cluster'].values          # 38 classes
y_global= query.obs['Global_Cluster'].values       # 5 broad classes

knn = KNeighborsClassifier(n_neighbors=15, metric='euclidean', n_jobs=-1)
knn.fit(X_gse, y_gse)
print(f"  KNN fitted on {len(X_gse):,} reference cells")

# ─────────────────────────────────────────────────────────────────────────────
# 4. Annotate immune cells in training set
# ─────────────────────────────────────────────────────────────────────────────
print("\n[4/6] Predicting sub-clusters for training set immune cells ...")

# Only annotate immune cells (those with matching types in GSE146771)
immune_types = {'T cells', 'B cells', 'NK cells', 'Myeloids'}
is_immune    = ref.obs['unified_cell_type'].isin(immune_types)
n_immune     = is_immune.sum()
print(f"  Immune cells to annotate: {n_immune:,}")

X_train_immune = ref.obsm['X_scvi'][is_immune.values]

# Predict
y_pred_sub  = knn.predict(X_train_immune)
y_prob      = knn.predict_proba(X_train_immune)
confidence  = y_prob.max(axis=1)

# Derive Global_Cluster label from sub-cluster prefix
prefix_to_global = {
    'hT': 'CD4/CD8 T cell',
    'hB': 'B cell',
    'hM': 'Myeloid cell',
    'hI': 'ILC/NK cell',
}
y_pred_global = np.array([
    prefix_to_global.get(sc[:2], 'Unknown') for sc in y_pred_sub
])

print(f"\n  Confidence stats:")
print(f"    mean={confidence.mean():.3f}  median={np.median(confidence):.3f}  "
      f"min={confidence.min():.3f}  max={confidence.max():.3f}")
print(f"    >= 0.5: {(confidence >= 0.5).sum():,}  "
      f"< 0.5:  {(confidence < 0.5).sum():,}")

# Write back to adata
CONF_THRESHOLD = 0.4   # lower threshold for sub-cluster (more types = lower per-type prob)
ref.obs['subcluster']           = 'Unassigned'
ref.obs['subcluster_global']    = 'Unassigned'
ref.obs['subcluster_conf']      = np.nan

immune_idx = ref.obs.index[is_immune.values]
ref.obs.loc[immune_idx, 'subcluster']        = y_pred_sub
ref.obs.loc[immune_idx, 'subcluster_global'] = y_pred_global
ref.obs.loc[immune_idx, 'subcluster_conf']   = confidence

high_conf_mask = confidence >= CONF_THRESHOLD
low_conf_idx   = immune_idx[~high_conf_mask]
ref.obs.loc[low_conf_idx, 'subcluster']        = 'Low-confidence'
ref.obs.loc[low_conf_idx, 'subcluster_global'] = 'Low-confidence'

print(f"\n  Sub-cluster distribution (training immune cells):")
sc_counts = ref.obs.loc[immune_idx, 'subcluster'].value_counts()
for name, n in sc_counts.items():
    pct = n / n_immune * 100
    print(f"    {name:35s}: {n:>6,}  ({pct:.1f}%)")

# ─────────────────────────────────────────────────────────────────────────────
# 5. Self-consistency check
# ─────────────────────────────────────────────────────────────────────────────
print("\n[5/6] Self-consistency check ...")

# Broad label vs sub-cluster prefix agreement
broad_map = {
    'T cells':  'hT',
    'B cells':  'hB',
    'Myeloids': 'hM',
    'NK cells': 'hI',
}
immune_df = ref.obs.loc[immune_idx].copy()
immune_df['expected_prefix'] = immune_df['unified_cell_type'].map(broad_map)
immune_df['actual_prefix']   = immune_df['subcluster'].str[:2]
immune_df['consistent']      = (
    (immune_df['subcluster'] == 'Low-confidence') |
    (immune_df['actual_prefix'] == immune_df['expected_prefix'])
)
consistent_high = (
    (immune_df['subcluster'] != 'Low-confidence') &
    (immune_df['actual_prefix'] == immune_df['expected_prefix'])
)

print(f"  High-conf predictions:          {high_conf_mask.sum():,}")
print(f"  Consistent with broad label:    {consistent_high.sum():,} "
      f"({consistent_high.sum()/high_conf_mask.sum()*100:.1f}% of high-conf)")
print()
print("  Consistency by broad type:")
for btype in ['T cells', 'B cells', 'Myeloids', 'NK cells']:
    sub  = immune_df[immune_df['unified_cell_type'] == btype]
    hc   = sub[sub['subcluster'] != 'Low-confidence']
    if len(hc) == 0:
        continue
    cons = (hc['actual_prefix'] == broad_map[btype]).sum()
    print(f"    {btype:12s}: {cons:,}/{len(hc):,} = {cons/len(hc)*100:.1f}% consistent")

# Cross-tab: our broad type vs predicted sub-cluster global
print()
print("  Cross-tab (broad label vs sub-cluster Global_Cluster):")
ct = pd.crosstab(
    immune_df['unified_cell_type'],
    immune_df['subcluster_global']
)
print(ct.to_string())

# ─────────────────────────────────────────────────────────────────────────────
# 6. Visualise & save
# ─────────────────────────────────────────────────────────────────────────────
print("\n[6/6] Generating plots & saving ...")

# --- Plot 1: UMAP coloured by sub-cluster (immune only shown) ----------------
fig, axes = plt.subplots(1, 3, figsize=(24, 8))

# All cells, sub-cluster (non-immune = grey)
plot_col = ref.obs['subcluster'].copy().astype(str)
plot_col[~is_immune.values] = 'Non-immune'
ref.obs['_plot_sub'] = plot_col

sc.pl.umap(ref, color='_plot_sub', ax=axes[0], show=False,
           title='Sub-cluster (GSE146771 labels)', frameon=False,
           legend_loc='right margin', legend_fontsize=6,
           palette=sns.color_palette('tab20', n_colors=40))
del ref.obs['_plot_sub']

# Sub-cluster confidence
sc.pl.umap(ref, color='subcluster_conf', ax=axes[1], show=False,
           title='Sub-cluster Confidence', frameon=False,
           cmap='RdYlGn', vmin=0, vmax=1)

# Broad immune type for comparison
sc.pl.umap(ref, color='unified_cell_type', ax=axes[2], show=False,
           title='Broad Cell Type', frameon=False,
           legend_loc='right margin', legend_fontsize=8)

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/umap_subclusters.png', dpi=200, bbox_inches='tight')
plt.close()
print(f"  Saved: {OUT_DIR}/umap_subclusters.png")

# --- Plot 2: Sub-cluster bar chart by broad type ----------------------------
fig, ax = plt.subplots(figsize=(14, 8))
high_conf_df = immune_df[immune_df['subcluster'] != 'Low-confidence'].copy()
ct2 = pd.crosstab(high_conf_df['subcluster'],
                  high_conf_df['unified_cell_type'])
# Sort sub-clusters by prefix then name
ct2 = ct2.sort_index()
ct2.plot(kind='barh', stacked=True, ax=ax,
         color=['#E64B35', '#4DBBD5', '#00A087', '#3C5488'],
         edgecolor='white', linewidth=0.3)
ax.set_xlabel('Cells', fontsize=11)
ax.set_title('Sub-cluster composition\n(high-confidence predictions)', fontsize=13, fontweight='bold')
ax.legend(title='Broad type', bbox_to_anchor=(1.01, 1), fontsize=9)
plt.tight_layout()
plt.savefig(f'{OUT_DIR}/subcluster_composition.png', dpi=200, bbox_inches='tight')
plt.close()
print(f"  Saved: {OUT_DIR}/subcluster_composition.png")

# --- Plot 3: Consistency heatmap -------------------------------------------
fig, ax = plt.subplots(figsize=(10, 5))
ct_norm = ct.div(ct.sum(axis=1), axis=0)
sns.heatmap(ct_norm, annot=ct, fmt='d', cmap='Blues', ax=ax,
            linewidths=0.5, cbar_kws={'shrink': 0.7})
ax.set_title('Broad label vs Sub-cluster Global\n(row normalised)',
             fontsize=12, fontweight='bold')
ax.set_xlabel('Predicted Sub-cluster Global', fontsize=10)
ax.set_ylabel('Our Broad Label', fontsize=10)
plt.tight_layout()
plt.savefig(f'{OUT_DIR}/consistency_heatmap.png', dpi=200, bbox_inches='tight')
plt.close()
print(f"  Saved: {OUT_DIR}/consistency_heatmap.png")

# --- Save annotated data ----------------------------------------------------
ref.write_h5ad('processed_data/merged_annotated.h5ad', compression='gzip')
print(f"  Updated: processed_data/merged_annotated.h5ad")

# Save sub-cluster table
sc_table = ref.obs[['batch', 'disease', 'unified_cell_type',
                     'subcluster', 'subcluster_global', 'subcluster_conf']].copy()
sc_table.to_csv(f'{OUT_DIR}/subcluster_assignments.csv')
print(f"  Saved: {OUT_DIR}/subcluster_assignments.csv")

print("\n" + "=" * 70)
print(f"COMPLETE — {datetime.now().strftime('%H:%M:%S')}")
n_assigned = (ref.obs['subcluster'] != 'Unassigned').sum()
n_highconf = (ref.obs.loc[is_immune.values[ref.obs.index.isin(immune_idx)],
              'subcluster'] != 'Low-confidence').sum() if False else high_conf_mask.sum()
print(f"  Immune cells annotated: {n_immune:,}")
print(f"  High-confidence (≥{CONF_THRESHOLD}): {high_conf_mask.sum():,} "
      f"({high_conf_mask.sum()/n_immune*100:.1f}%)")
print("=" * 70)
