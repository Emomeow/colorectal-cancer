#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 2: Train scVI on the 3000-HVG merged dataset.

Input : processed_data/merged_hvg3k.h5ad  (195,090 cells x ~3,000 genes)
Output: scvi_model_hvg3k/
        scvi_results_hvg3k/merged_scvi_hvg3k.h5ad
"""

import scanpy as sc
import scvi
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import os
import sys
import io

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

print("=" * 70)
print(f"[{datetime.now().strftime('%H:%M:%S')}] scVI Training on HVG-3k subset")
print("=" * 70)

DATA_PATH  = 'processed_data/merged_hvg3k.h5ad'
MODEL_DIR  = 'scvi_model_hvg3k'
RESULTS_DIR = 'scvi_results_hvg3k'

N_LATENT  = 30
N_LAYERS  = 2
MAX_EPOCHS = 400
LEARNING_RATE = 1e-3
BATCH_SIZE = 256
EARLY_STOPPING_PATIENCE = 50

os.makedirs(MODEL_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── Step 1: Load data ─────────────────────────────────────────────────────────
print(f"\n[1] Loading {DATA_PATH} ...")
adata = sc.read_h5ad(DATA_PATH)
print(f"    {adata.n_obs:,} cells x {adata.n_vars:,} genes")
print(f"    Batches: {adata.obs['batch'].nunique()}")
print(f"    Diseases: {list(adata.obs['disease'].unique())}")
print(f"    layers: {list(adata.layers.keys())}")

# ── Step 2: Setup scVI ────────────────────────────────────────────────────────
print(f"\n[2] Setting up scVI ...")
scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts',                          # raw integer counts
    batch_key='batch',
    categorical_covariate_keys=['disease'],
)
print(f"    Setup complete (using layer='counts')")

# ── Step 3: Initialize model ──────────────────────────────────────────────────
print(f"\n[3] Initializing model ...")
model = scvi.model.SCVI(
    adata,
    n_latent=N_LATENT,
    n_layers=N_LAYERS,
    gene_likelihood='nb',
)
print(f"    {adata.n_vars:,} genes -> {N_LAYERS} layers -> {N_LATENT} latent dims")

# ── Step 4: Train ─────────────────────────────────────────────────────────────
print(f"\n[4] Training ...")
print(f"    lr={LEARNING_RATE}  batch_size={BATCH_SIZE}  max_epochs={MAX_EPOCHS}  "
      f"early_stop_patience={EARLY_STOPPING_PATIENCE}")

start_time = datetime.now()
print(f"    Started: {start_time.strftime('%H:%M:%S')}")

model.train(
    max_epochs=MAX_EPOCHS,
    batch_size=BATCH_SIZE,
    early_stopping=True,
    early_stopping_patience=EARLY_STOPPING_PATIENCE,
    plan_kwargs={
        'lr': LEARNING_RATE,
        'reduce_lr_on_plateau': True,
        'lr_scheduler_metric': 'elbo_validation',
        'lr_patience': 8,
        'lr_factor': 0.6,
    },
)

end_time = datetime.now()
duration = (end_time - start_time).total_seconds()
print(f"    Finished: {end_time.strftime('%H:%M:%S')}  ({duration/60:.1f} min)")

# ── Step 5: Training history ──────────────────────────────────────────────────
print(f"\n[5] Training history ...")
train_history = model.history
try:
    if isinstance(train_history, pd.DataFrame):
        hist = train_history
    elif hasattr(train_history, '__getitem__'):
        # dict-like history (scvi >= 0.20)
        keys = list(train_history.keys())
        hist = pd.DataFrame({k: train_history[k]['elbo_train'].values
                              if hasattr(train_history[k], 'values')
                              else train_history[k]
                              for k in keys if 'elbo' in k})
    else:
        hist = None

    if hist is not None and 'elbo_train' in hist.columns:
        print(f"    Epochs trained: {len(hist)}")
        print(f"    Final train ELBO:      {hist['elbo_train'].iloc[-1]:.2f}")
        if 'elbo_validation' in hist.columns:
            print(f"    Final val   ELBO:      {hist['elbo_validation'].iloc[-1]:.2f}")
    else:
        hist = None
        print(f"    History format: {type(train_history)}")
except Exception as e:
    print(f"    Warning: {e}")
    hist = None

# ── Step 6: Save model ────────────────────────────────────────────────────────
print(f"\n[6] Saving model to {MODEL_DIR}/ ...")
model.save(MODEL_DIR, overwrite=True)
model_mb = sum(os.path.getsize(os.path.join(MODEL_DIR, f))
               for f in os.listdir(MODEL_DIR)
               if os.path.isfile(os.path.join(MODEL_DIR, f))) / 1024**2
print(f"    Model size: {model_mb:.1f} MB")

# ── Step 7: Latent embeddings ─────────────────────────────────────────────────
print(f"\n[7] Generating latent embeddings ...")
latent = model.get_latent_representation()
adata.obsm['X_scvi'] = latent
print(f"    Embeddings: {latent.shape}")

# ── Step 8: Neighbor graph + UMAP ────────────────────────────────────────────
print(f"\n[8] Neighbor graph + UMAP ...")
sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=15)
sc.tl.umap(adata, min_dist=0.3)
print(f"    UMAP done")

# ── Step 9: Save h5ad ─────────────────────────────────────────────────────────
out_h5ad = os.path.join(RESULTS_DIR, 'merged_scvi_hvg3k.h5ad')
print(f"\n[9] Saving integrated data to {out_h5ad} ...")
adata.write_h5ad(out_h5ad, compression='gzip')
print(f"    {os.path.getsize(out_h5ad)/1024**2:.1f} MB")

# ── Step 10: Training history plot ───────────────────────────────────────────
print(f"\n[10] Plotting training history ...")
if hist is not None:
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype']  = 42
    plt.rcParams['font.family']  = 'Arial'

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    epochs = range(1, len(hist) + 1)

    if 'elbo_train' in hist.columns:
        axes[0].plot(epochs, hist['elbo_train'].values, label='Train', lw=2)
        if 'elbo_validation' in hist.columns:
            axes[0].plot(epochs, hist['elbo_validation'].values, label='Val', lw=2)
        axes[0].set(xlabel='Epoch', ylabel='ELBO', title='ELBO')
        axes[0].legend(); axes[0].grid(alpha=0.3)

    if 'reconstruction_loss_train' in hist.columns:
        axes[1].plot(epochs, hist['reconstruction_loss_train'].values, label='Train', lw=2)
        if 'reconstruction_loss_validation' in hist.columns:
            axes[1].plot(epochs, hist['reconstruction_loss_validation'].values, label='Val', lw=2)
        axes[1].set(xlabel='Epoch', ylabel='Reconstruction Loss', title='Reconstruction Loss')
        axes[1].legend(); axes[1].grid(alpha=0.3)

    plt.tight_layout()
    for ext in ['png', 'pdf']:
        p = os.path.join(RESULTS_DIR, f'training_history.{ext}')
        plt.savefig(p, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"    Saved training_history.png/.pdf")
else:
    print(f"    Skipped (history unavailable)")

# ── Step 11: UMAP plots ───────────────────────────────────────────────────────
print(f"\n[11] UMAP plots ...")
if 'X_umap' in adata.obsm:
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype']  = 42
    plt.rcParams['font.family']  = 'Arial'

    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    sc.pl.umap(adata, color='batch',             ax=axes[0, 0], show=False,
               title='Batch', frameon=False)
    sc.pl.umap(adata, color='disease',           ax=axes[0, 1], show=False,
               title='Disease', frameon=False,
               palette={'CRC': '#E64B35', 'UC': '#4DBBD5'})
    sc.pl.umap(adata, color='unified_cell_type', ax=axes[1, 0], show=False,
               title='Cell Type', frameon=False, legend_loc='right margin')
    sc.pl.umap(adata, color='n_genes_by_counts', ax=axes[1, 1], show=False,
               title='Gene Count', frameon=False, cmap='viridis')
    plt.tight_layout()
    for ext in ['png', 'pdf']:
        p = os.path.join(RESULTS_DIR, f'umap_overview.{ext}')
        plt.savefig(p, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"    Saved umap_overview.png/.pdf")

# ── Summary ───────────────────────────────────────────────────────────────────
print(f"\n{'='*70}")
print(f"DONE")
print(f"  Cells:    {adata.n_obs:,}")
print(f"  Genes:    {adata.n_vars:,} (HVG-3k)")
print(f"  Latent:   {N_LATENT}D")
print(f"  Duration: {duration/60:.1f} min")
print(f"  Model:    {MODEL_DIR}/")
print(f"  Data:     {out_h5ad}")
print(f"[{datetime.now().strftime('%H:%M:%S')}] Complete.")
