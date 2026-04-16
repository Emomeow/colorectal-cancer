#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Assemble Figure 1 from individual panel PNGs."""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from PIL import Image
import os

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42
plt.rcParams['font.family'] = 'Arial'

PANELS = {
    'AB': 'figures/phase1/fig1_panelAB.png',
    'C':  'figures/phase1/fig1_panelC.png',
    'D':  'figures/phase1/fig1_panelD.png',
    'EF': 'figures/phase1/fig1_panelEF.png',
}

# Check all files exist
for k, p in PANELS.items():
    if not os.path.exists(p):
        raise FileNotFoundError(f"Missing panel: {p}")

# Load as RGB (drop alpha if present) — required for PDF embedding
imgs = {k: np.array(Image.open(p).convert('RGB')) for k, p in PANELS.items()}

# ── Layout: 2 rows
#   Row 1: AB (wider) | C
#   Row 2: D          | EF (wider)
# Use gridspec width_ratios to balance
fig = plt.figure(figsize=(20, 13))
gs = gridspec.GridSpec(
    2, 2,
    figure=fig,
    hspace=0.06,
    wspace=0.04,
    left=0.01, right=0.99,
    top=0.97, bottom=0.01,
    width_ratios=[1.85, 1.0],   # AB/EF wider; C/D narrower
    height_ratios=[1.0, 1.0],
)

ax_ab = fig.add_subplot(gs[0, 0])
ax_c  = fig.add_subplot(gs[0, 1])
ax_d  = fig.add_subplot(gs[1, 0])
ax_ef = fig.add_subplot(gs[1, 1])

for ax, key in [(ax_ab, 'AB'), (ax_c, 'C'), (ax_d, 'D'), (ax_ef, 'EF')]:
    ax.imshow(imgs[key])
    ax.axis('off')

os.makedirs('figures/phase1', exist_ok=True)

# Save PNG
png_out = 'figures/phase1/figure1_landscape.png'
fig.savefig(png_out, dpi=200, bbox_inches='tight')
plt.close(fig)
size_mb = os.path.getsize(png_out) / 1024**2
print(f"Saved: {png_out}  ({size_mb:.1f} MB)")

# Convert PNG → PDF via PIL (more reliable than matplotlib PDF backend for raster content)
pdf_out = 'figures/phase1/figure1_landscape.pdf'
Image.open(png_out).convert('RGB').save(pdf_out, 'PDF', resolution=200)
size_mb = os.path.getsize(pdf_out) / 1024**2
print(f"Saved: {pdf_out}  ({size_mb:.1f} MB)")

print("Done.")
