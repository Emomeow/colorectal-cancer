#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXPERIMENTAL: scVI with HVG-union gene set
- Selects top-N HVGs per batch, takes union of genes present in >=MIN_DATASETS datasets
- Pads missing genes with zeros (marked as unmeasured)
- Results saved to scvi_hvgunion/ — original files NOT touched

Original model: scvi_model_optimized/  (9,009 intersection genes)
This experiment:  scvi_hvgunion/
"""
import scanpy as sc
import scvi
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys, io
from datetime import datetime

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

os.makedirs('scvi_hvgunion', exist_ok=True)

N_HVG_PER_BATCH = 3000   # HVGs to select per dataset
MIN_DATASETS    = 2       # gene must be measured in >= this many datasets
N_LATENT        = 30
N_LAYERS        = 2
MAX_EPOCHS      = 200
BATCH_SIZE      = 256

print("=" * 65)
print(f"HVG-union scVI experiment  |  {datetime.now().strftime('%H:%M:%S')}")
print(f"  N_HVG_PER_BATCH={N_HVG_PER_BATCH}  MIN_DATASETS={MIN_DATASETS}")
print("=" * 65)

# ─────────────────────────────────────────────────────────────────
# 1. Load individual datasets and compute per-batch HVGs
# ─────────────────────────────────────────────────────────────────
PATHS = {
    'GSE132257': ('processed_data/GSE132257.h5ad', 'CRC'),
    'GSE132465': ('processed_data/GSE132465.h5ad', 'CRC'),
    'GSE125527': ('processed_data/GSE125527.h5ad', 'UC'),
    'GSE150115': ('processed_data/GSE150115.h5ad', 'UC'),
    'cellxgene': ('processed_data/cellxgene_raw.h5ad', 'CRC'),
}

print("\n[1/6] Loading datasets and computing per-batch HVGs ...")
hvg_per_batch = {}
adatas_raw    = {}

for batch, (path, disease) in PATHS.items():
    adata = sc.read_h5ad(path)
    adata.obs['batch']   = batch
    adata.obs['disease'] = disease
    adatas_raw[batch]    = adata

    # HVG on normalized copy
    tmp = adata.copy()
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)
    sc.pp.highly_variable_genes(tmp, n_top_genes=N_HVG_PER_BATCH, flavor='seurat_v3',
                                 span=0.3)
    hvgs = tmp.var_names[tmp.var['highly_variable']].tolist()
    hvg_per_batch[batch] = set(hvgs)
    print(f"  {batch:12s}  cells={adata.n_obs:6,}  genes={adata.n_vars:6,}  HVGs={len(hvgs)}")

# ─────────────────────────────────────────────────────────────────
# 2. Build union gene set filtered by MIN_DATASETS
# ─────────────────────────────────────────────────────────────────
print(f"\n[2/6] Building HVG union (genes in >={MIN_DATASETS} datasets) ...")

# Count how many datasets each HVG gene appears in (measured, not just HVG)
all_measured = {b: set(adatas_raw[b].var_names) for b in PATHS}
all_hvg_union = set.union(*hvg_per_batch.values())

# For each gene in HVG union: how many datasets actually MEASURE it?
gene_counts = {}
for gene in all_hvg_union:
    n = sum(1 for b in PATHS if gene in all_measured[b])
    gene_counts[gene] = n

# Filter
union_genes = sorted([g for g, c in gene_counts.items() if c >= MIN_DATASETS])
print(f"  HVG union (all):              {len(all_hvg_union):,}")
for thresh in [2, 3, 4, 5]:
    n = sum(1 for c in gene_counts.values() if c >= thresh)
    print(f"  Present in >={thresh} datasets:       {n:,} genes")
print(f"\n  Using: {len(union_genes):,} genes  (MIN_DATASETS={MIN_DATASETS})")

# How many are NEW vs intersection?
original_genes = set(sc.read_h5ad('processed_data/merged_annotated.h5ad').var_names)
new_genes   = set(union_genes) - original_genes
lost_genes  = original_genes  - set(union_genes)
print(f"  vs original 9,009 intersection:")
print(f"    New genes added : {len(new_genes):,}")
print(f"    Original genes kept : {len(set(union_genes) & original_genes):,}")
print(f"    Original genes lost : {len(lost_genes):,}  (not HVG in any batch)")

# ─────────────────────────────────────────────────────────────────
# 3. Build merged dataset with zero-padding for missing genes
# ─────────────────────────────────────────────────────────────────
print(f"\n[3/6] Building merged dataset with {len(union_genes):,} genes ...")

def align_to_genes(X_src, src_genes, tgt_genes):
    """Reindex columns of X_src to tgt_genes, filling missing with zeros."""
    tgt_idx = {g: i for i, g in enumerate(tgt_genes)}
    n_cells, n_tgt = X_src.shape[0], len(tgt_genes)

    # Map source column j -> target column
    col_map = [(j, tgt_idx[g]) for j, g in enumerate(src_genes) if g in tgt_idx]
    if not col_map:
        return sp.csr_matrix((n_cells, n_tgt), dtype=np.float32)

    src_cols, tgt_cols = map(list, zip(*col_map))
    Xs = sp.csr_matrix(X_src)[:, src_cols].tocoo()

    result = sp.coo_matrix(
        (Xs.data, (Xs.row, np.array(tgt_cols)[Xs.col])),
        shape=(n_cells, n_tgt)
    ).tocsr()
    return result

adatas_aligned = []
for batch, (path, disease) in PATHS.items():
    adata = adatas_raw[batch]
    X_new = align_to_genes(adata.X, list(adata.var_names), union_genes)
    adata_new = sc.AnnData(
        X   = X_new,
        obs = adata.obs.copy(),
        var = pd.DataFrame(index=union_genes),
    )
    # Track which genes were actually measured in this batch
    measured_mask = np.array([g in all_measured[batch] for g in union_genes])
    adata_new.var[f'measured_{batch}'] = measured_mask
    adatas_aligned.append(adata_new)
    n_zero = (~measured_mask).sum()
    print(f"  {batch:12s}  zero-padded genes: {n_zero:,} / {len(union_genes):,}")

merged = sc.concat(adatas_aligned, join='outer', label='batch', keys=list(PATHS.keys()))
merged.obs_names_make_unique()
print(f"\n  Merged: {merged.n_obs:,} cells x {merged.n_vars:,} genes")

# ─────────────────────────────────────────────────────────────────
# 4. Train scVI
# ─────────────────────────────────────────────────────────────────
print(f"\n[4/6] Training scVI ...")
scvi.model.SCVI.setup_anndata(
    merged,
    batch_key='batch',
    categorical_covariate_keys=['disease'],
    layer=None,
)
model = scvi.model.SCVI(
    merged,
    n_latent=N_LATENT,
    n_layers=N_LAYERS,
    gene_likelihood='nb',
    dropout_rate=0.1,
)
print(f"  Architecture: {merged.n_vars:,} genes -> {N_LAYERS} layers -> {N_LATENT} latent dims")

model.train(
    max_epochs=MAX_EPOCHS,
    batch_size=BATCH_SIZE,
    early_stopping=True,
    early_stopping_patience=15,
    early_stopping_monitor='elbo_validation',
    train_size=0.9,
    plan_kwargs={'lr': 1e-3},
)
model.save('scvi_hvgunion/model', overwrite=True)
print("  Model saved to scvi_hvgunion/model/")

# ─────────────────────────────────────────────────────────────────
# 5. Get latent + UMAP
# ─────────────────────────────────────────────────────────────────
print(f"\n[5/6] Computing latent representations and UMAP ...")
merged.obsm['X_scvi'] = model.get_latent_representation()
sc.pp.neighbors(merged, use_rep='X_scvi', n_neighbors=30)
sc.tl.umap(merged, min_dist=0.3)
merged.write_h5ad('scvi_hvgunion/merged_hvgunion.h5ad')
print("  Saved: scvi_hvgunion/merged_hvgunion.h5ad")

# ─────────────────────────────────────────────────────────────────
# 6. Quick comparison plots: UMAP by batch and by cell type
# ─────────────────────────────────────────────────────────────────
print(f"\n[6/6] Generating comparison UMAP ...")

# Load original for comparison
orig = sc.read_h5ad('processed_data/merged_annotated.h5ad')

BATCH_COLORS = {
    'GSE132257': '#E64B35', 'GSE132465': '#4DBBD5',
    'GSE125527': '#00A087', 'GSE150115': '#3C5488',
    'cellxgene': '#F39B7F',
}
TYPE_COLORS = {
    'T cells': '#E64B35', 'NK cells': '#3C5488', 'B cells': '#4DBBD5',
    'Myeloids': '#00A087', 'Mast cells': '#F39B7F', 'Epithelial': '#91D1C2',
    'Stromal': '#DC0000', 'Malignant': '#B09C85', 'Unknown': '#AAAAAA',
}

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

def umap_scatter(ax, adata, color_col, color_map, title, s=1.5, alpha=0.4):
    for label, color in color_map.items():
        mask = adata.obs[color_col] == label
        if mask.sum() == 0:
            continue
        xy = adata.obsm['X_umap'][mask.values]
        ax.scatter(xy[:, 0], xy[:, 1], c=color, s=s, alpha=alpha, label=label, linewidths=0)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.axis('off')
    ax.legend(markerscale=4, fontsize=7, loc='lower left', framealpha=0.7)

# Top row: original (9,009 genes)
umap_scatter(axes[0, 0], orig, 'batch', BATCH_COLORS,
             f'Original  ({len(original_genes):,} genes)\nby dataset')
umap_scatter(axes[0, 1], orig, 'cell_type_grouped', TYPE_COLORS,
             f'Original  ({len(original_genes):,} genes)\nby cell type')

# Bottom row: HVG-union
umap_scatter(axes[1, 0], merged, 'batch', BATCH_COLORS,
             f'HVG-union ({len(union_genes):,} genes, min={MIN_DATASETS} datasets)\nby dataset')

# For cell type: merged doesn't have cell_type_grouped yet — use unified from obs if available
ct_col = 'cell_type_grouped' if 'cell_type_grouped' in merged.obs.columns else 'batch'
umap_scatter(axes[1, 1], merged, ct_col, TYPE_COLORS if ct_col == 'cell_type_grouped' else BATCH_COLORS,
             f'HVG-union ({len(union_genes):,} genes)\nby {ct_col}')

plt.suptitle('HVG-union vs Original: UMAP comparison', fontsize=13, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig('scvi_hvgunion/umap_comparison.png', dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved: scvi_hvgunion/umap_comparison.png")

# Training loss curve
train_hist = model.history['elbo_train']
val_hist   = model.history.get('elbo_validation', None)
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(train_hist, label='train ELBO')
if val_hist is not None:
    ax.plot(val_hist.values, label='val ELBO')
ax.set_xlabel('Epoch'); ax.set_ylabel('ELBO'); ax.legend()
ax.set_title(f'HVG-union scVI training  ({len(union_genes):,} genes)')
plt.tight_layout()
plt.savefig('scvi_hvgunion/training_curve.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved: scvi_hvgunion/training_curve.png")

print(f"\n{'='*65}")
print(f"Done  |  {datetime.now().strftime('%H:%M:%S')}")
print(f"  Gene set: {len(union_genes):,} (vs original 9,009)")
print(f"  Results: scvi_hvgunion/")
print(f"  Original results untouched.")
print(f"{'='*65}")
