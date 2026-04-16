#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Optimized scVI training with better parameters and faster convergence
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

# Set encoding for Windows
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

print("=" * 80)
print(f"[{datetime.now().strftime('%H:%M:%S')}] Optimized scVI Training")
print("=" * 80)

# ============================================================================
# OPTIMIZED CONFIGURATION
# ============================================================================

DATA_PATH = 'processed_data/merged_5datasets_unified.h5ad'
MODEL_DIR = 'scvi_model_optimized'
RESULTS_DIR = 'scvi_results_optimized'

# Model architecture
N_LATENT = 30           # Latent space dimensions (20-40 recommended)
N_LAYERS = 2            # Hidden layers (1-2 recommended)

# Training parameters - OPTIMIZED
MAX_EPOCHS = 100        # Reduced from 200 (since it plateaus early)
LEARNING_RATE = 1e-3    # Increased from default 4e-4
BATCH_SIZE = 256        # Increased from default 128 for speed
EARLY_STOPPING = True
EARLY_STOPPING_PATIENCE = 15  # Reduced from 45

print(f"\nOptimized Configuration:")
print(f"  Latent dimensions: {N_LATENT}")
print(f"  Hidden layers: {N_LAYERS}")
print(f"  Max epochs: {MAX_EPOCHS}")
print(f"  Learning rate: {LEARNING_RATE}")
print(f"  Batch size: {BATCH_SIZE}")
print(f"  Early stopping patience: {EARLY_STOPPING_PATIENCE}")

# Create output directories
os.makedirs(MODEL_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# ============================================================================
# STEP 1: Load Data
# ============================================================================
print(f"\n{'='*80}")
print(f"[Step 1] Loading Data")
print(f"{'='*80}")

try:
    adata = sc.read_h5ad(DATA_PATH)
    print(f"Data loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"Batches: {len(adata.obs['batch'].unique())}")
    print(f"Diseases: {', '.join(adata.obs['disease'].unique())}")
except Exception as e:
    print(f"Error loading data: {e}")
    sys.exit(1)

# ============================================================================
# STEP 2: Setup scVI
# ============================================================================
print(f"\n{'='*80}")
print(f"[Step 2] Setting up scVI")
print(f"{'='*80}")

try:
    scvi.model.SCVI.setup_anndata(
        adata,
        batch_key='batch',
        categorical_covariate_keys=['disease'],
        layer=None
    )
    print(f"scVI setup complete")
except Exception as e:
    print(f"Error setting up scVI: {e}")
    sys.exit(1)

# ============================================================================
# STEP 3: Initialize Model
# ============================================================================
print(f"\n{'='*80}")
print(f"[Step 3] Initializing Model")
print(f"{'='*80}")

try:
    model = scvi.model.SCVI(
        adata,
        n_latent=N_LATENT,
        n_layers=N_LAYERS,
        gene_likelihood='nb'
    )
    print(f"Model initialized")
    print(f"  Architecture: {adata.n_vars:,} genes -> {N_LAYERS} layers -> {N_LATENT} latent dims")
except Exception as e:
    print(f"Error initializing model: {e}")
    sys.exit(1)

# ============================================================================
# STEP 4: Train Model (OPTIMIZED)
# ============================================================================
print(f"\n{'='*80}")
print(f"[Step 4] Training Model (Optimized)")
print(f"{'='*80}")

print(f"\nTraining with optimized parameters:")
print(f"  Learning rate: {LEARNING_RATE}")
print(f"  Batch size: {BATCH_SIZE}")
print(f"  Max epochs: {MAX_EPOCHS}")
print(f"  Early stopping: {EARLY_STOPPING} (patience={EARLY_STOPPING_PATIENCE})")

start_time = datetime.now()
print(f"\nTraining started: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

try:
    model.train(
        max_epochs=MAX_EPOCHS,
        batch_size=BATCH_SIZE,
        early_stopping=EARLY_STOPPING,
        early_stopping_patience=EARLY_STOPPING_PATIENCE,
        plan_kwargs={
            'lr': LEARNING_RATE,
            'reduce_lr_on_plateau': True,
            'lr_scheduler_metric': 'elbo_validation',
            'lr_patience': 8,
            'lr_factor': 0.6
        }
    )

    end_time = datetime.now()
    duration = (end_time - start_time).total_seconds()

    print(f"\nTraining complete!")
    print(f"  End time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Duration: {duration/60:.1f} minutes")

except Exception as e:
    print(f"\nTraining failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# ============================================================================
# STEP 5: Get Training History (FIXED)
# ============================================================================
print(f"\n{'='*80}")
print(f"[Step 5] Getting Training History")
print(f"{'='*80}")

try:
    train_history = model.history

    if isinstance(train_history, pd.DataFrame):
        if 'elbo_train' in train_history.columns:
            # FIX: Use .iloc[-1] instead of [-1] for DataFrame
            final_elbo = train_history['elbo_train'].iloc[-1]
            print(f"Final training ELBO: {final_elbo:.2f}")

            if 'elbo_validation' in train_history.columns:
                final_val_elbo = train_history['elbo_validation'].iloc[-1]
                print(f"Final validation ELBO: {final_val_elbo:.2f}")

            # Calculate per-cell loss
            per_cell_loss = final_elbo / adata.n_obs
            print(f"Per-cell average loss: {per_cell_loss:.4f}")

            print(f"\nTraining epochs: {len(train_history)}")
    else:
        print(f"Training history format: {type(train_history)}")

except Exception as e:
    print(f"Warning: Could not get training history: {e}")
    train_history = None

# ============================================================================
# STEP 6: Save Model
# ============================================================================
print(f"\n{'='*80}")
print(f"[Step 6] Saving Model")
print(f"{'='*80}")

try:
    model.save(MODEL_DIR, overwrite=True)

    saved_files = os.listdir(MODEL_DIR)
    model_size = sum(os.path.getsize(os.path.join(MODEL_DIR, f))
                     for f in saved_files if os.path.isfile(os.path.join(MODEL_DIR, f)))

    print(f"Model saved: {MODEL_DIR}")
    print(f"  Size: {model_size / (1024**2):.1f} MB")

except Exception as e:
    print(f"Error saving model: {e}")
    import traceback
    traceback.print_exc()

# ============================================================================
# STEP 7: Generate Embeddings
# ============================================================================
print(f"\n{'='*80}")
print(f"[Step 7] Generating Latent Embeddings")
print(f"{'='*80}")

try:
    print(f"Generating embeddings (this may take a few minutes)...")
    latent = model.get_latent_representation()

    print(f"Embeddings generated: {latent.shape[0]:,} cells x {latent.shape[1]} dims")

    adata.obsm['X_scvi'] = latent
    print(f"Added to adata.obsm['X_scvi']")

except Exception as e:
    print(f"Error generating embeddings: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# ============================================================================
# STEP 8: Save Integrated Data
# ============================================================================
print(f"\n{'='*80}")
print(f"[Step 8] Saving Integrated Data")
print(f"{'='*80}")

try:
    output_path = os.path.join(RESULTS_DIR, 'merged_scvi_integrated.h5ad')
    adata.write_h5ad(output_path, compression='gzip')

    file_size = os.path.getsize(output_path) / (1024**2)
    print(f"Integrated data saved: {output_path}")
    print(f"  Size: {file_size:.1f} MB")

except Exception as e:
    print(f"Error saving integrated data: {e}")
    import traceback
    traceback.print_exc()

# ============================================================================
# STEP 9: Plot Training History
# ============================================================================
print(f"\n{'='*80}")
print(f"[Step 9] Plotting Training History")
print(f"{'='*80}")

try:
    if train_history is not None and isinstance(train_history, pd.DataFrame):
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))

        epochs = range(1, len(train_history) + 1)

        # Plot 1: ELBO
        if 'elbo_train' in train_history.columns:
            axes[0].plot(epochs, train_history['elbo_train'].values,
                        label='Training', linewidth=2)
            if 'elbo_validation' in train_history.columns:
                axes[0].plot(epochs, train_history['elbo_validation'].values,
                           label='Validation', linewidth=2)

            axes[0].set_xlabel('Epoch', fontsize=12)
            axes[0].set_ylabel('ELBO', fontsize=12)
            axes[0].set_title('Training History - ELBO', fontsize=14, fontweight='bold')
            axes[0].legend()
            axes[0].grid(True, alpha=0.3)

        # Plot 2: Reconstruction Loss
        if 'reconstruction_loss_train' in train_history.columns:
            axes[1].plot(epochs, train_history['reconstruction_loss_train'].values,
                        label='Training', linewidth=2)
            if 'reconstruction_loss_validation' in train_history.columns:
                axes[1].plot(epochs, train_history['reconstruction_loss_validation'].values,
                           label='Validation', linewidth=2)

            axes[1].set_xlabel('Epoch', fontsize=12)
            axes[1].set_ylabel('Reconstruction Loss', fontsize=12)
            axes[1].set_title('Training History - Reconstruction Loss',
                            fontsize=14, fontweight='bold')
            axes[1].legend()
            axes[1].grid(True, alpha=0.3)

        plt.tight_layout()

        plot_path = os.path.join(RESULTS_DIR, 'training_history.png')
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"Training history plot saved: {plot_path}")
    else:
        print(f"Skipping training history plot (history not available)")

except Exception as e:
    print(f"Warning: Could not plot training history: {e}")

# ============================================================================
# STEP 10: Compute UMAP
# ============================================================================
print(f"\n{'='*80}")
print(f"[Step 10] Computing UMAP")
print(f"{'='*80}")

try:
    print(f"Computing neighbor graph...")
    sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=15)
    print(f"Neighbor graph computed")

    print(f"Computing UMAP...")
    sc.tl.umap(adata, min_dist=0.3)
    print(f"UMAP computed")

    # Save with UMAP
    output_path = os.path.join(RESULTS_DIR, 'merged_scvi_integrated.h5ad')
    adata.write_h5ad(output_path, compression='gzip')
    print(f"UMAP added and saved")

except Exception as e:
    print(f"Error computing UMAP: {e}")
    import traceback
    traceback.print_exc()

# ============================================================================
# STEP 11: Generate UMAP Visualizations
# ============================================================================
print(f"\n{'='*80}")
print(f"[Step 11] Generating UMAP Visualizations")
print(f"{'='*80}")

try:
    if 'X_umap' in adata.obsm:
        fig, axes = plt.subplots(2, 2, figsize=(16, 14))

        # Plot 1: Batch
        sc.pl.umap(adata, color='batch', ax=axes[0, 0], show=False,
                   title='scVI Integration - Batch', frameon=False)

        # Plot 2: Disease
        sc.pl.umap(adata, color='disease', ax=axes[0, 1], show=False,
                   title='scVI Integration - Disease', frameon=False,
                   palette={'CRC': '#E64B35', 'UC': '#4DBBD5'})

        # Plot 3: Cell Type
        sc.pl.umap(adata, color='unified_cell_type', ax=axes[1, 0], show=False,
                   title='scVI Integration - Cell Type', frameon=False,
                   legend_loc='right margin')

        # Plot 4: Gene Count
        sc.pl.umap(adata, color='n_genes_by_counts', ax=axes[1, 1], show=False,
                   title='scVI Integration - Gene Count', frameon=False,
                   cmap='viridis')

        plt.tight_layout()

        plot_path = os.path.join(RESULTS_DIR, 'umap_overview.png')
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"UMAP visualization saved: {plot_path}")

except Exception as e:
    print(f"Error generating UMAP visualization: {e}")
    import traceback
    traceback.print_exc()

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print(f"\n{'='*80}")
print(f"TRAINING COMPLETE")
print(f"{'='*80}")

print(f"\nGenerated Files:")
print(f"  1. Model: {MODEL_DIR}/")
print(f"  2. Integrated data: {RESULTS_DIR}/merged_scvi_integrated.h5ad")
print(f"  3. Training history: {RESULTS_DIR}/training_history.png")
print(f"  4. UMAP visualization: {RESULTS_DIR}/umap_overview.png")

print(f"\nTraining Summary:")
print(f"  Duration: {duration/60:.1f} minutes")
print(f"  Epochs: {len(train_history) if train_history is not None else 'N/A'}")
print(f"  Latent dimensions: {N_LATENT}")
print(f"  Learning rate: {LEARNING_RATE}")
print(f"  Batch size: {BATCH_SIZE}")

print(f"\nNext Steps:")
print(f"  1. Evaluate batch mixing (kBET, LISI)")
print(f"  2. Cell type annotation for Unknown cells")
print(f"  3. Differential analysis (CRC vs UC)")
print(f"  4. Marker gene analysis")

print(f"\n{'='*80}")
print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"{'='*80}")
