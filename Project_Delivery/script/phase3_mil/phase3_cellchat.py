#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 3 — Cell Communication Analysis
OLR1+ TAM (SPP1+B) <-> FAP+ Myofibroblast ligand-receptor interactions
Dataset: GSE144735 (27,414 cells, Border/Normal/Tumor)
Method: LIANA (CellPhoneDB + Magnitude consensus)
"""
import sys, io, os, warnings
import pandas as pd
import numpy as np
import anndata
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
warnings.filterwarnings('ignore')
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import liana
from liana.method import cellphonedb, rank_aggregate

plt.rcParams['font.family']  = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42

OUT_CSV = 'phase3_results/cellchat/'
OUT_FIG = 'Project_Delivery/results/Figure3_MIL_MMR_Prediction/figures/'
os.makedirs(OUT_CSV, exist_ok=True)
os.makedirs(OUT_FIG, exist_ok=True)

# ── Load GSE144735 ─────────────────────────────────────────────────────
print("Loading GSE144735...")
adata = anndata.read_h5ad('processed_data/GSE144735.h5ad')
print(f"  Total cells: {adata.n_obs}, genes: {adata.n_vars}")

# Normalize if not already
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ── Analysis 1: All Border cells (context for interaction) ─────────────
print("\n--- Analysis 1: Border-region communication ---")
border = adata[adata.obs['Class'] == 'Border'].copy()
print(f"  Border cells: {border.n_obs}")
print(border.obs['Cell_subtype'].value_counts()[
    ['SPP1+B', 'Myofibroblasts', 'CD8+ T cells', 'CD4+ T cells',
     'Regulatory T cells', 'Anti-inflammatory', 'Pro-inflammatory']
])

# Run CellPhoneDB on Border region
# Rename Cell_subtype for liana (needs 'groupby' column)
border.obs['cell_type'] = border.obs['Cell_subtype'].astype(str)

print("\nRunning LIANA rank_aggregate on Border cells...")
try:
    rank_aggregate(border, groupby='cell_type',
                   expr_prop=0.05, min_cells=5,
                   verbose=True, n_perms=100,
                   use_raw=False)
    liana_res_border = border.uns['liana_res'].copy()
    print(f"  Total interactions found: {len(liana_res_border)}")
except Exception as e:
    print(f"rank_aggregate failed: {e}, trying cellphonedb...")
    cellphonedb(border, groupby='cell_type',
                expr_prop=0.05, min_cells=5,
                n_perms=100, verbose=True, use_raw=False)
    liana_res_border = border.uns['liana_res'].copy()
    print(f"  Total interactions found: {len(liana_res_border)}")

liana_res_border.to_csv(f'{OUT_CSV}liana_border_all.csv', index=False)

# ── Analysis 2: SPP1+B -> Myofibroblast (key story pair) ──────────────
print("\n--- Filtering: SPP1+B <-> Myofibroblasts ---")

# Forward: SPP1+B sends to Myofibroblasts
spp1_to_fibro = liana_res_border[
    (liana_res_border['source'] == 'SPP1+B') &
    (liana_res_border['target'] == 'Myofibroblasts')
].copy()

# Reverse: Myofibroblasts sends to SPP1+B
fibro_to_spp1 = liana_res_border[
    (liana_res_border['source'] == 'Myofibroblasts') &
    (liana_res_border['target'] == 'SPP1+B')
].copy()

print(f"  SPP1+B -> Myofibroblasts: {len(spp1_to_fibro)} pairs")
print(f"  Myofibroblasts -> SPP1+B: {len(fibro_to_spp1)} pairs")

# ── Sort by magnitude score ────────────────────────────────────────────
# liana rank_aggregate produces 'magnitude_rank' (lower = more significant)
# or cellphonedb produces 'lr_means' and 'cellphone_pvals'
score_col = None
for c in ['magnitude_rank', 'lr_means', 'lrs_to_filter', 'mean']:
    if c in spp1_to_fibro.columns:
        score_col = c
        break

print(f"  Scoring column: {score_col}")
print(f"  Columns: {spp1_to_fibro.columns.tolist()}")

if score_col == 'magnitude_rank':
    spp1_to_fibro_top = spp1_to_fibro.sort_values(score_col).head(20)
    fibro_to_spp1_top = fibro_to_spp1.sort_values(score_col).head(15)
else:
    # fallback: sort by lr_means descending
    alt = [c for c in spp1_to_fibro.columns if 'mean' in c.lower() or 'score' in c.lower()]
    if alt:
        spp1_to_fibro_top = spp1_to_fibro.sort_values(alt[0], ascending=False).head(20)
        fibro_to_spp1_top = fibro_to_spp1.sort_values(alt[0], ascending=False).head(15)
    else:
        spp1_to_fibro_top = spp1_to_fibro.head(20)
        fibro_to_spp1_top = fibro_to_spp1.head(15)

spp1_to_fibro_top.to_csv(f'{OUT_CSV}spp1b_to_myofibro_top.csv', index=False)
fibro_to_spp1_top.to_csv(f'{OUT_CSV}myofibro_to_spp1b_top.csv', index=False)

print("\nTop SPP1+B -> Myofibroblasts pairs:")
print(spp1_to_fibro_top[['ligand_complex', 'receptor_complex'] +
      [score_col] if score_col else
      spp1_to_fibro_top.iloc[:, :5]].to_string())

# ── Analysis 3: Compare Border vs Tumor communication ──────────────────
print("\n--- Analysis 3: Tumor-region communication for comparison ---")
tumor = adata[adata.obs['Class'] == 'Tumor'].copy()
tumor.obs['cell_type'] = tumor.obs['Cell_subtype'].astype(str)

# Check if SPP1+B and Myofibroblasts exist in Tumor
tumor_cts = tumor.obs['Cell_subtype'].value_counts()
print(f"  SPP1+B in Tumor: {tumor_cts.get('SPP1+B', 0)}")
print(f"  Myofibroblasts in Tumor: {tumor_cts.get('Myofibroblasts', 0)}")

if tumor_cts.get('SPP1+B', 0) >= 5 and tumor_cts.get('Myofibroblasts', 0) >= 5:
    try:
        rank_aggregate(tumor, groupby='cell_type',
                       expr_prop=0.05, min_cells=5,
                       verbose=False, n_perms=100, use_raw=False)
        liana_res_tumor = tumor.uns['liana_res'].copy()
        spp1_tumor = liana_res_tumor[
            (liana_res_tumor['source'] == 'SPP1+B') &
            (liana_res_tumor['target'] == 'Myofibroblasts')
        ]
        spp1_tumor.to_csv(f'{OUT_CSV}spp1b_to_myofibro_tumor.csv', index=False)
        print(f"  SPP1+B->Myofibro interactions in Tumor: {len(spp1_tumor)}")
    except Exception as e:
        print(f"  Tumor analysis error: {e}")

# ── Plot Figure: Bubble plot of top L-R pairs ──────────────────────────
print("\nGenerating communication figures...")

def plot_lr_bubble(df_fwd, df_rev, score_col, out_prefix):
    """Bubble plot: x=lr_pair, y=direction, size=score, color=ligand_type"""
    # Determine value column for sizing
    if score_col == 'magnitude_rank':
        # Convert rank to score (lower rank = higher importance)
        df_fwd = df_fwd.copy()
        df_rev = df_rev.copy()
        max_rank = max(df_fwd[score_col].max() if len(df_fwd) else 1,
                       df_rev[score_col].max() if len(df_rev) else 1)
        df_fwd['plot_score'] = 1 - df_fwd[score_col] / (max_rank + 1)
        df_rev['plot_score'] = 1 - df_rev[score_col] / (max_rank + 1)
    else:
        # Use mean expression score
        val_col = [c for c in df_fwd.columns if 'mean' in c.lower()]
        val_col = val_col[0] if val_col else df_fwd.columns[-1]
        df_fwd = df_fwd.copy()
        df_rev = df_rev.copy()
        df_fwd['plot_score'] = df_fwd[val_col]
        df_rev['plot_score'] = df_rev[val_col]

    # pval column
    pval_col = None
    for c in df_fwd.columns:
        if 'pval' in c.lower() or 'p_val' in c.lower() or 'pvalue' in c.lower():
            pval_col = c
            break

    # Story pairs to highlight
    STORY_PAIRS = {
        'SPP1_CD44', 'SPP1_ITGAV', 'SPP1_ITGB1', 'SPP1_ITGA4',
        'MMP9_LRP1', 'MMP9_ITGB2', 'TREM2_APOE',
        'CCL20_CCR6', 'INHBA_ACVR1', 'FAP_',
        'SPP1', 'MMP9', 'TREM2', 'OLR1',
        'CD44', 'FAP', 'ACTA2', 'COL1A1',
    }

    fig, axes = plt.subplots(1, 2, figsize=(18, 7),
                              gridspec_kw={'width_ratios': [1.3, 1.0]})
    fig.suptitle(
        'OLR1+ TAM (SPP1+B) <-> FAP+ Myofibroblast\nCell Communication (Border Region, GSE144735)',
        fontsize=12, fontweight='bold', y=1.02
    )

    for ax, df, title, direction_color in [
        (axes[0], df_fwd, 'A  OLR1+ TAM → Myofibroblast\n(TAM as Sender)', '#E64B35'),
        (axes[1], df_rev, 'B  Myofibroblast → OLR1+ TAM\n(Fibro as Sender)', '#2CA02C'),
    ]:
        if len(df) == 0:
            ax.text(0.5, 0.5, 'No interactions\nfound', ha='center', va='center',
                    fontsize=11, color='gray', transform=ax.transAxes)
            ax.set_title(title, fontsize=10, fontweight='bold', loc='left')
            continue

        df = df.head(20 if ax == axes[0] else 15).copy().reset_index(drop=True)
        pairs = df['ligand_complex'].astype(str) + ' → ' + df['receptor_complex'].astype(str)

        scores = df['plot_score'].values
        scores = np.clip(scores, 0, 1)

        # Color: highlight story pairs
        colors = []
        for _, row in df.iterrows():
            lig = str(row['ligand_complex'])
            rec = str(row['receptor_complex'])
            is_story = any(s in lig or s in rec for s in STORY_PAIRS)
            colors.append('#E64B35' if is_story else '#AAAAAA')

        y_pos = np.arange(len(df))
        sizes = 20 + scores * 300

        for y, fc, col, s in zip(y_pos, scores, colors, sizes):
            ax.plot([0, fc], [y, y], color=col, lw=1.5, alpha=0.7,
                    solid_capstyle='round')
            ax.scatter(fc, y, s=s, color=col, alpha=0.85,
                       edgecolors='black', linewidths=0.6, zorder=5)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(pairs.tolist(), fontsize=8, fontstyle='italic')

        # Bold story pairs
        for lbl, (_, row) in zip(ax.get_yticklabels(), df.iterrows()):
            lig = str(row['ligand_complex'])
            rec = str(row['receptor_complex'])
            if any(s in lig or s in rec for s in STORY_PAIRS):
                lbl.set_fontweight('bold')
                lbl.set_color('#CC3300')

        ax.axvline(0, color='#888888', lw=0.8, ls='--')
        ax.set_xlabel('Interaction Score', fontsize=9)
        ax.set_title(title, fontsize=10, fontweight='bold', loc='left')
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='x', color='#EEEEEE', lw=0.5)

    story_patch = mpatches.Patch(color='#E64B35', label='Story L-R pair (SPP1/MMP9/TREM2→CAF)')
    other_patch = mpatches.Patch(color='#AAAAAA', label='Other interactions')
    axes[0].legend(handles=[story_patch, other_patch], fontsize=8,
                   loc='lower right', framealpha=0.85)

    fig.text(0.5, -0.01,
             'LIANA rank_aggregate (CellPhoneDB + OmniPath methods) | Border region only (n=9,424 cells) | '
             'Size ∝ interaction score | Bold red = OLR1+ TAM core story pairs',
             ha='center', fontsize=7, color='#555555')

    plt.tight_layout()
    for ext in ['png', 'pdf']:
        fig.savefig(f'{OUT_FIG}{out_prefix}.{ext}',
                    dpi=200 if ext == 'png' else None,
                    bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {OUT_FIG}{out_prefix}.png / .pdf")

plot_lr_bubble(spp1_to_fibro_top, fibro_to_spp1_top, score_col,
               'cellchat_olr1tam_myofibro')

# ── Summary table for report ───────────────────────────────────────────
print("\n--- Summary tables ---")
print("\nTop 10 SPP1+B -> Myofibroblasts:")
cols_show = ['ligand_complex', 'receptor_complex']
if score_col:
    cols_show.append(score_col)
if 'cellphone_pvals' in spp1_to_fibro_top.columns:
    cols_show.append('cellphone_pvals')
print(spp1_to_fibro_top[cols_show].head(10).to_string(index=False))

print("\nTop 10 Myofibroblasts -> SPP1+B:")
print(fibro_to_spp1_top[cols_show].head(10).to_string(index=False))

print("\nDone.")
