#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 1 – Panel C: Cell-type composition per batch (CellAssign annotations)"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42
plt.rcParams['font.family'] = 'Arial'
os.makedirs('figures/phase1', exist_ok=True)

adata = sc.read_h5ad('scvi_results_hvg3k/merged_scvi_hvg3k.h5ad')

# CellAssign cell types (10 types)
CT_ORDER = ['CD4 T', 'CD8 T', 'Treg', 'NK',
            'B cell', 'Plasma',
            'Myeloid', 'Mast',
            'Stromal', 'Epithelial']

CT_COLORS = {
    'CD4 T':     '#4DBBD5',
    'CD8 T':     '#00A087',
    'Treg':      '#3C5488',
    'NK':        '#F39B7F',
    'B cell':    '#F7DC6F',
    'Plasma':    '#8491B4',
    'Myeloid':   '#91D1C2',
    'Mast':      '#E64B35',
    'Stromal':   '#B09C85',
    'Epithelial':'#BBCC33',
}

# ── Composition for all 5 batches ─────────────────────────────────────────────
comp = (adata.obs
        .groupby(['batch', 'cellass_type'])
        .size()
        .unstack(fill_value=0))

# Ensure all cell types present
for ct in CT_ORDER:
    if ct not in comp.columns:
        comp[ct] = 0
comp = comp[CT_ORDER]

BATCH_ORDER = ['GSE132257', 'GSE132465', 'cellxgene', 'GSE125527', 'GSE150115']
comp = comp.reindex(BATCH_ORDER, fill_value=0)
comp_pct = comp.div(comp.sum(axis=1), axis=0) * 100

print("Composition (%) per batch:")
print(comp_pct.round(1).to_string())

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9.5, 4.8))
fig.subplots_adjust(left=0.08, right=0.74, top=0.90, bottom=0.16)

x   = np.arange(5)
W   = 0.72
bottoms = np.zeros(5)

for ct in CT_ORDER:
    vals = comp_pct[ct].values
    ax.bar(x, vals, bottom=bottoms, color=CT_COLORS[ct],
           width=W, linewidth=0.3, edgecolor='white', label=ct)
    bottoms += vals

# Cell counts below bars
for xi, batch in enumerate(BATCH_ORDER):
    n = comp.loc[batch].sum()
    ax.text(xi, -4, f'n={n:,}', ha='center', va='top',
            fontsize=6.5, color='#444')

# Disease labels above bars
DISEASE_LABEL = {
    'GSE132257': ('CRC', '#E64B35'),
    'GSE132465': ('CRC', '#E64B35'),
    'cellxgene': ('CRC', '#E64B35'),
    'GSE125527': ('UC',  '#4DBBD5'),
    'GSE150115': ('UC',  '#4DBBD5'),
}
for xi, batch in enumerate(BATCH_ORDER):
    label, color = DISEASE_LABEL[batch]
    ax.text(xi, 102, label, ha='center', fontsize=8.5,
            color=color, fontweight='bold')

# Separator between CRC and UC
ax.axvline(2.6, color='#999', lw=0.8, ls=':', alpha=0.7)

XTICKLABELS = ['GSE132257', 'GSE132465', 'CXG', 'GSE125527', 'GSE150115']
ax.set_xticks(x)
ax.set_xticklabels(XTICKLABELS, fontsize=8)
ax.set_ylabel('Proportion (%)', fontsize=9)
ax.set_xlim(-0.6, 4.6)
ax.set_ylim(-10, 110)
ax.tick_params(labelsize=7.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.text(-0.07, 1.03, 'D', transform=ax.transAxes,
        fontsize=14, fontweight='bold', va='bottom')
ax.set_title('Cell Type Composition per Dataset (CellAssign)', fontsize=10, pad=4)

leg = [Patch(color=CT_COLORS[ct], label=ct) for ct in CT_ORDER]
ax.legend(handles=leg, fontsize=7.5, frameon=False, ncol=1,
          loc='center left', bbox_to_anchor=(1.01, 0.5),
          handlelength=0.9, handleheight=0.9, labelspacing=0.4)

for ext in ['png', 'pdf']:
    fig.savefig(f'figures/phase1/fig1_panelC.{ext}',
                dpi=200, bbox_inches='tight')
plt.close(fig)
print("Saved figures/phase1/fig1_panelC.png/.pdf")
