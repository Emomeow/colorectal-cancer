#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rebuild merged_annotated.h5ad from scratch with correct priority order:

Priority 1 (highest): unified_cell_subtype already in merged h5ad
  → Rescues NK cells (1,111) and T/B subtypes from GSE132465 original paper

Priority 2: GSE125527 metadata file (T/B/NK/M_DC)
  → Covers UC immune cells with original paper labels

Priority 3: KNN label transfer via scVI embeddings
  → Fills remaining Unknown cells

Also adds:
  - annotation_source column tracking every cell's label origin
  - Fine-grained subtype labels where available
"""

import scanpy as sc
import scvi
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.neighbors import KNeighborsClassifier
import os, sys, io
from datetime import datetime

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

OUT_DIR = 'annotation_results'
os.makedirs(OUT_DIR, exist_ok=True)

t0 = datetime.now()
print("=" * 70)
print("Rebuilding merged_annotated.h5ad — full annotation pipeline")
print("=" * 70)

# ─────────────────────────────────────────────────────────────────────────────
# 1. Load original (intact) merged data + scVI results
# ─────────────────────────────────────────────────────────────────────────────
print("\n[1/7] Loading original merged data ...")
adata = sc.read_h5ad('processed_data/merged_5datasets_unified.h5ad')

# Load scVI embeddings and UMAP from the optimised run (these are intact)
scvi_res = sc.read_h5ad('scvi_results_optimized/merged_scvi_integrated.h5ad')
adata.obsm['X_scvi'] = scvi_res.obsm['X_scvi']
adata.obsm['X_umap'] = scvi_res.obsm['X_umap']

print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")
print(f"  X_scvi: {adata.obsm['X_scvi'].shape}")
print(f"  X_umap: {adata.obsm['X_umap'].shape}")

# Working copy of unified_cell_type
adata.obs['unified_cell_type'] = adata.obs['unified_cell_type'].astype(str)
# Track annotation source
src = pd.Series('original', index=adata.obs.index, dtype=str)

# ─────────────────────────────────────────────────────────────────────────────
# 2. Priority 1 — Rescue annotations from unified_cell_subtype (GSE132465)
# ─────────────────────────────────────────────────────────────────────────────
print("\n[2/7] Priority 1 — Rescuing subtype annotations (GSE132465 original) ...")

# Subtype → broad type mapping for subtypes that change the parent label
subtype_upgrade = {
    # NK cells: currently under "T cells", should be "NK cells"
    'NK cells':               'NK cells',
    # T cell subtypes — keep unified_cell_type as T cells, store subtype
    'CD4+ T cells':           'T cells',
    'CD8+ T cells':           'T cells',
    'Regulatory T cells':     'T cells',
    'T follicular helper cells': 'T cells',
    'T helper 17 cells':      'T cells',
    'gamma delta T cells':    'T cells',
    # B cell subtypes
    'CD19+CD20+ B':           'B cells',
    'IgA+ Plasma':            'B cells',
    'IgG+ Plasma':            'B cells',
    # Myeloid subtypes
    'SPP1+':                  'Myeloids',
    'Pro-inflammatory':       'Myeloids',
    'cDC':                    'Myeloids',
    'Proliferating':          'Myeloids',
}

sub_col = adata.obs['unified_cell_subtype'].astype(str)
for subtype, broad in subtype_upgrade.items():
    mask = sub_col == subtype
    if mask.sum() == 0:
        continue
    old = adata.obs.loc[mask, 'unified_cell_type'].value_counts().to_dict()
    adata.obs.loc[mask, 'unified_cell_type'] = broad
    src[mask] = f'original_subtype ({subtype})'
    print(f"  {subtype:40s} -> {broad:12s}: {mask.sum():,} cells  (was: {old})")

# ─────────────────────────────────────────────────────────────────────────────
# 3. Priority 2 — GSE125527 metadata recovery (T/B/NK/M_DC)
# ─────────────────────────────────────────────────────────────────────────────
print("\n[3/7] Priority 2 — GSE125527 metadata recovery ...")

meta   = pd.read_csv('data/GSE125527_cell_metadata.csv.gz')
id_map = pd.read_csv('data/GSE125527_oldPatientId-newPatientId.csv.gz')
old2new = dict(zip(id_map['old_id'], id_map['new_id']))

celltype_map = {
    'T':       'T cells',
    'B':       'B cells',
    'M/DC':    'Myeloids',
    'NK':      'NK cells',
    'unknown': None,   # skip — leave as Unknown
}
meta['unified_type'] = meta['celltype'].map(celltype_map)
meta = meta[meta['unified_type'].notna()]

meta_composite = (meta['patient_assignment'] + '|' + meta['cell_id']).values
meta_type_by_c = dict(zip(meta_composite, meta['unified_type']))

gse_mask = adata.obs['batch'] == 'GSE125527'
gse_cells = adata.obs[gse_mask].copy()
stripped  = gse_cells.index.str.replace('_GSE125527', '', regex=False)
parts     = stripped.str.split('_', n=2)
gse_cells['new_patient'] = parts.str[0].map(old2new)
gse_cells['raw_barcode'] = parts.str[2]
gse_cells['composite']   = gse_cells['new_patient'].fillna('') + '|' + gse_cells['raw_barcode']

new_types    = gse_cells['composite'].map(meta_type_by_c)
matched_mask = new_types.notna()
n_meta = matched_mask.sum()

recovered_idx = gse_cells.index[matched_mask]
# Only overwrite cells that are still Unknown (don't overwrite Priority 1)
still_unknown = adata.obs.loc[recovered_idx, 'unified_cell_type'] == 'Unknown'
final_idx = recovered_idx[still_unknown]
adata.obs.loc[final_idx, 'unified_cell_type'] = new_types[matched_mask][still_unknown].values
src[final_idx] = 'GSE125527_metadata'

print(f"  Matched in metadata: {n_meta:,}")
print(f"  Written (was Unknown): {still_unknown.sum():,}")
cts = new_types[matched_mask][still_unknown].value_counts()
for ct, n in cts.items():
    print(f"    {ct}: {n:,}")

# ─────────────────────────────────────────────────────────────────────────────
# 4. Priority 3 — KNN label transfer for remaining Unknown cells
# ─────────────────────────────────────────────────────────────────────────────
print("\n[4/7] Priority 3 — KNN label transfer for remaining Unknown ...")

is_unknown = adata.obs['unified_cell_type'].isin(['Unknown', 'Unspecified'])
is_labeled = ~is_unknown
n_unknown  = is_unknown.sum()
print(f"  Labeled (reference): {is_labeled.sum():,}")
print(f"  Unknown (query):     {n_unknown:,}")

X_ref   = adata.obsm['X_scvi'][is_labeled.values]
y_ref   = adata.obs.loc[is_labeled, 'unified_cell_type'].values
X_query = adata.obsm['X_scvi'][is_unknown.values]

knn = KNeighborsClassifier(n_neighbors=15, metric='euclidean', n_jobs=-1)
knn.fit(X_ref, y_ref)
y_pred      = knn.predict(X_query)
confidence  = knn.predict_proba(X_query).max(axis=1)

CONF = 0.5
high_conf = confidence >= CONF
unknown_idx = adata.obs.index[is_unknown.values]
adata.obs.loc[unknown_idx[high_conf], 'unified_cell_type'] = y_pred[high_conf]
src[unknown_idx[high_conf]] = 'KNN_transfer'
src[unknown_idx[~high_conf]] = 'KNN_low_confidence'

# Store KNN metadata
adata.obs['knn_predicted_type'] = 'N/A'
adata.obs['knn_confidence']     = np.nan
adata.obs.loc[unknown_idx, 'knn_predicted_type'] = y_pred
adata.obs.loc[unknown_idx, 'knn_confidence']     = confidence

unknown_final = adata.obs['unified_cell_type'].isin(['Unknown', 'Unspecified']).sum()
print(f"  High-conf KNN (>={CONF}): {high_conf.sum():,}")
print(f"  Still Unknown:            {unknown_final:,} ({unknown_final/adata.n_obs*100:.1f}%)")
knn_cts = pd.Series(y_pred[high_conf]).value_counts()
for ct, n in knn_cts.items():
    print(f"    {ct}: {n:,}")

# ─────────────────────────────────────────────────────────────────────────────
# 5. Finalise + annotation_source column
# ─────────────────────────────────────────────────────────────────────────────
print("\n[5/7] Finalising columns ...")

adata.obs['unified_cell_type']  = adata.obs['unified_cell_type'].astype('category')
adata.obs['annotation_source']  = src.astype('category')

# NK cell reliability flag
nk_mask = adata.obs['unified_cell_type'] == 'NK cells'
nk_src  = adata.obs.loc[nk_mask, 'annotation_source']
print(f"\n  NK cells total: {nk_mask.sum():,}")
print("  NK cells by source:")
for s, n in nk_src.value_counts().items():
    flag = '✓ reliable' if 'subtype' in str(s) or 'metadata' in str(s) else '⚠ inferred'
    print(f"    {s}: {n:,}  {flag}")

print("\n  Final unified_cell_type distribution:")
for ct, n in adata.obs['unified_cell_type'].value_counts().items():
    print(f"    {ct:45s}: {n:,}")

# ─────────────────────────────────────────────────────────────────────────────
# 6. Save
# ─────────────────────────────────────────────────────────────────────────────
print("\n[6/7] Saving ...")
out = 'processed_data/merged_annotated.h5ad'
adata.write_h5ad(out, compression='gzip')
size = os.path.getsize(out) / (1024**2)
print(f"  Saved: {out}  ({size:.1f} MB)")

# ─────────────────────────────────────────────────────────────────────────────
# 7. Summary plots
# ─────────────────────────────────────────────────────────────────────────────
print("\n[7/7] Summary plots ...")

fig, axes = plt.subplots(1, 3, figsize=(24, 8))

sc.pl.umap(adata, color='unified_cell_type', ax=axes[0], show=False,
           title='Cell Type (Final)', frameon=False,
           legend_loc='right margin', legend_fontsize=7)

sc.pl.umap(adata, color='annotation_source', ax=axes[1], show=False,
           title='Annotation Source', frameon=False,
           legend_loc='right margin', legend_fontsize=7)

sc.pl.umap(adata, color='batch', ax=axes[2], show=False,
           title='Batch', frameon=False,
           legend_loc='right margin', legend_fontsize=8)

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/umap_final_annotation.png', dpi=200, bbox_inches='tight')
plt.close()
print(f"  Saved: {OUT_DIR}/umap_final_annotation.png")

# Annotation source summary table
src_ct = pd.crosstab(adata.obs['annotation_source'],
                     adata.obs['unified_cell_type'])
src_ct.to_csv(f'{OUT_DIR}/annotation_source_by_celltype.csv')
print(f"  Saved: {OUT_DIR}/annotation_source_by_celltype.csv")

elapsed = (datetime.now() - t0).total_seconds()
print(f"\n{'='*70}")
print(f"COMPLETE in {elapsed/60:.1f} min — {datetime.now().strftime('%H:%M:%S')}")
print(f"  NK cells from original paper (GSE132465): "
      f"{(nk_src.str.contains('subtype')).sum():,}")
print(f"  NK cells from GSE125527 metadata:         "
      f"{(nk_src.str.contains('metadata')).sum():,}")
print(f"  NK cells from KNN (low reliability):      "
      f"{(nk_src.str.contains('KNN')).sum():,}")
print(f"{'='*70}")
